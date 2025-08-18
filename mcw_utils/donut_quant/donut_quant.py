from transformers import AutoImageProcessor, ResNetModel
from datasets import load_dataset
import torch

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


def cos_sim(a, b):
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))


def get_sample(row):
    sample, _ = row["index"].split("_")
    return sample


def get_drug(row):
    _, drug = row["index"].split("_")
    return drug


def combine_replicates(mean_cos_sims, indices_df, file_df):
    exp_time = file_df.time[1:]
    indices_df_out = indices_df.copy()
    for t in exp_time:
        indices_df_out[t] = -1.0
    for t in exp_time:
        indices_df_out[f"{t}_stdev"] = -1.0
    for cell_line, drug in indices_df.index:
        stack = []
        for i, j in indices_df.loc[cell_line].loc[drug]:
            stack.append(np.array(mean_cos_sims[j][i]))
        mean_cos_sim = np.mean(stack, axis=0)
        std_cos_sim = np.std(stack, axis=0)
        for t, m in zip(exp_time, mean_cos_sim):
            indices_df_out.loc[(cell_line, drug), t] = m
        for t, s in zip(exp_time, std_cos_sim):
            indices_df_out.loc[(cell_line, drug), f"{t}_stdev"] = s
    return indices_df_out


def compute_cos_sim(time_encodings, indices_df):
    vehicle_label = None
    for index in indices_df.index:
        sample, drug = index
        if drug == "vehicle":
            vehicle_label = "vehicle"
        elif drug == 0:
            vehicle_label = 0
    if vehicle_label is None:
        raise ValueError(
            "No vehicle label found in drug.csv, must be either 'vehicle' or 0"
        )
    vehicle_df = indices_df.swaplevel(0, 1).loc[vehicle_label, :].copy()
    average_zero_tensor_matrix = {j: {} for j in range(12)}
    mean_cos_sims = {j: {} for j in range(12)}
    # std_cos_sims = {j: {} for j in range(12)}

    for sample in vehicle_df.index:
        sample_df = indices_df.loc[sample, :]

        for replicate in sample_df.columns:
            stack = []
            # print(sample_df[replicate])
            for i, j in sample_df[replicate]:
                for replicate_rotation_encoding in time_encodings[0][j][i]:
                    # for rotation in encoding:
                    stack.append(replicate_rotation_encoding)
                stacked_tensors = torch.stack(stack)
                average_zero_tensor_matrix[j][i] = torch.mean(stacked_tensors, dim=0)
    for sample in vehicle_df.index:
        sample_df = indices_df.loc[sample, :]
        for replicate in sample_df.columns:
            stack = []
            # print(sample_df[replicate])
            for i, j in sample_df[replicate]:
                mean_cos_sim_time = []
                # std_cos_sim_time = []
                for time in range(len(time_encodings)):
                    cos_sims = []
                    for replicate_rotation_encoding in time_encodings[time][j][i]:
                        # for rotation in encoding:
                        cos_sims.append(
                            cos_sim(
                                average_zero_tensor_matrix[j][i],
                                replicate_rotation_encoding,
                            )
                        )
                    mean_cos_sim_time.append(np.mean(cos_sims))
                    # std_cos_sim_time.append(np.std(cos_sims))
                mean_cos_sims[j][i] = mean_cos_sim_time
                # std_cos_sims[j][i] = std_cos_sim_time
    for sample in vehicle_df.index:
        sample_df = indices_df.loc[sample, :]
        for replicate in sample_df.columns:
            stack = []
            # print(sample_df[replicate])
            for i, j in sample_df[replicate]:
                stack.append(mean_cos_sims[j][i])

    return mean_cos_sims


def crop_and_encode(
    img_file, l=53, u=51, w=325, step=447, show_plot=False, rotations=4
):
    """
    smol = l=112, u=112, w=200, step=448
    large =w=325, l=53,u=51, step=447
    Measure the vector representations of images at each time point relative to vehicle.

    This crops a 96 well plate image at maximum LICOR resolution
    """
    encodings = []
    axs = None
    if show_plot:
        _, axs = plt.subplots(8, 12)

    # img_file = folder_path / f"{st}.png"  # Replace with your image URL or path
    image = Image.open(img_file).convert("L").convert("RGB")
    image_processor = AutoImageProcessor.from_pretrained("microsoft/resnet-50")
    model = ResNetModel.from_pretrained("microsoft/resnet-50")
    for j in range(12):
        encoding = []
        for i in range(8):
            cropped = image.crop(
                (l + j * step, u + i * step, l + j * step + w, u + i * step + w)
            )
            if show_plot:
                axs[i, j].imshow(cropped)
                axs[i, j].xaxis.set_visible(False)
                axs[i, j].yaxis.set_visible(False)
            rotation_encoding = []
            for i in range(rotations):
                cropped = cropped.transpose(Image.ROTATE_90)
                inputs = image_processor(cropped, return_tensors="pt")
                with torch.no_grad():
                    outputs = model(**inputs)
                    rotation_encoding.append(outputs.pooler_output.squeeze())
            encoding.append(rotation_encoding)

        encodings.append(encoding)
    return encodings


def read_base_path(base_path):
    file_df = pd.read_csv(base_path / "files.csv")

    file_names = []
    for folder in file_df.folder:
        folder_path = base_path / folder
        for f in folder_path.iterdir():
            if f.is_dir():
                tif_found = False
                for inner_files in f.iterdir():
                    if (
                        inner_files.is_file()
                        and inner_files.suffix == ".TIF"
                        and tif_found == False
                    ):
                        file_names.append(inner_files)
                        tif_found = True
                if tif_found == False:
                    raise ValueError(
                        f"No TIF file found in folder {f}. Please check the folder structure."
                    )
    file_df["tif_path"] = file_names
    sample_df = pd.read_csv(base_path / "samples.csv", header=None)
    drug_df = pd.read_csv(base_path / "drug.csv", header=None)
    sample_df.fillna(-1.0, inplace=True)
    drug_df.fillna(-1.0, inplace=True)
    sample_drug_combos = []
    for i in range(len(sample_df)):
        for s, d in zip(sample_df.iloc[i], drug_df.iloc[i]):
            if s != -1.0 and d != -1.0:
                sample_drug_combos.append(f"{s}_{d}")
    sample_drug_combos = list(set(sample_drug_combos))

    indices = {}
    for sample_drug_combo in sample_drug_combos:
        index = []
        for i in range(len(sample_df)):
            for j, sd in enumerate(zip(sample_df.iloc[i], drug_df.iloc[i])):
                s, d = sd
                if f"{s}_{d}" == sample_drug_combo:
                    index.append([i, j])
        indices[sample_drug_combo] = index

    indices_df = pd.DataFrame(indices)
    indices_df = indices_df.T.reset_index()
    indices_df["sample"] = indices_df.apply(get_sample, axis=1)
    indices_df["drug"] = indices_df.apply(get_drug, axis=1)
    indices_df.drop(columns=["index"], inplace=True)
    indices_df.sort_values(["sample", "drug"], inplace=True)
    indices_df.set_index(["sample", "drug"], inplace=True)
    for column in indices_df.columns:
        indices_df[f"replicate_{column}"] = indices_df[column]
        indices_df.drop(columns=[column], inplace=True)
    return file_df, indices_df


def main(base_path, show_plot=False):

    print("Reading inputs...")
    file_df, indices_df = read_base_path(base_path)
    print("Cropping and encoding wells...")
    time_encodings = []
    for i in range(len(file_df)):
        print(f"Processing {file_df.iloc[i].tif_path} at time {file_df.iloc[i].time}")
        encodings = crop_and_encode(file_df.iloc[i].tif_path, show_plot=show_plot)
        time_encodings.append(encodings)
    print("Encodings complete, computing cosine similarities...")
    mean_cos_sims = compute_cos_sim(time_encodings, indices_df)
    print("Combining replicates...")
    out_df = combine_replicates(mean_cos_sims, indices_df, file_df)
    out_df.to_csv(base_path / "resnet50_cosine_similarity.csv")
    print(f"Results saved to: {base_path / 'resnet50_cosine_similarity.csv'}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Measure donuts in file path.")
    parser.add_argument(
        "folder_path",
        type=str,
        help="Path to the folder containing cropped donut images.",
    )
    args = parser.parse_args()
    folder_path = Path(args.folder_path)
    main(folder_path, show_plot=False)
