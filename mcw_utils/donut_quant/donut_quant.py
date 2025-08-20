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
    return float(drug)


def compute_cos_sim(encodings, indices_df):
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

    sample_mean = {}
    sample_std = {}
    for sample in vehicle_df.index:
        vehicle_sample_df = vehicle_df.loc[sample, :]
        stack = []
        if isinstance(vehicle_sample_df, pd.Series):
            # If the sample is a Series, convert it to a DataFrame
            vehicle_sample_df = vehicle_sample_df.to_frame().T

        for replicate in vehicle_sample_df.columns:
            for i, j in vehicle_sample_df[replicate]:
                for replicate_rotation_encoding in encodings[j][i]:
                    # for rotation in encoding:
                    stack.append(replicate_rotation_encoding)
        stacked_tensors = torch.stack(stack)
        average_zero_tensor = torch.mean(stacked_tensors, dim=0)

        mean_cos_sims = {}
        std_cos_sims = {}
        sample_df = indices_df.loc[sample, :]
        for drug in sample_df.index:
            cos_sims = []

            drug_df = sample_df.loc[drug]
            stack = []
            if isinstance(drug_df, pd.Series):
                # If the sample is a Series, convert it to a DataFrame
                drug_df = drug_df.to_frame().T

            cos_sims = []
            for replicate in drug_df.columns:
                for i, j in drug_df[replicate]:
                    for replicate_rotation_encoding in encodings[j][i]:
                        # for rotation in encoding:
                        cos_sims.append(
                            cos_sim(average_zero_tensor, replicate_rotation_encoding)
                        )
            mean_cos_sims[drug] = np.mean(cos_sims)
            std_cos_sims[drug] = np.std(cos_sims)
        sample_mean[sample] = mean_cos_sims
        sample_std[sample] = std_cos_sims
    sample_mean_df = pd.DataFrame(sample_mean)
    sample_std_df = pd.DataFrame(sample_std)
    out_df = sample_mean_df.merge(
        sample_std_df,
        left_index=True,
        right_index=True,
        how="inner",
        suffixes=["_mean", "_std"],
    )
    return out_df


def crop_and_encode(
    img_file,
    l=53,
    u=51,
    w=325,
    step=447,
    plot_show=False,
    rotations=4,
    rows=8,
    columns=12,
):
    """
    smol = l=112, u=112, w=200, step=448
    large =w=325, l=53,u=51, step=447
    Measure the vector representations of images at each time point relative to vehicle.

    This crops a 96 well plate image at maximum LICOR resolution
    """
    encodings = []
    _, axs = plt.subplots(rows, columns)

    # img_file = folder_path / f"{st}.png"  # Replace with your image URL or path
    image = Image.open(img_file).convert("RGB")
    image_processor = AutoImageProcessor.from_pretrained("microsoft/resnet-50")
    model = ResNetModel.from_pretrained("microsoft/resnet-50")
    for j in range(columns):
        encoding = []
        for i in range(rows):
            cropped = image.crop(
                (l + j * step, u + i * step, l + j * step + w, u + i * step + w)
            )
            if plot_show:
                axs[i, j].imshow(cropped)
                axs[i, j].xaxis.set_visible(False)
                axs[i, j].yaxis.set_visible(False)
            rotation_encoding = []
            for i in range(rotations):
                cropped = cropped.transpose(Image.ROTATE_90)
                inputs = image_processor(cropped, return_tensors="pt", use_fast=True)
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
    rows = sample_df.shape[0]
    columns = sample_df.shape[1]
    sample_df.fillna("", inplace=True)
    drug_df.fillna("", inplace=True)
    sample_drug_combos = []
    for i in range(len(sample_df)):
        for s, d in zip(sample_df.iloc[i], drug_df.iloc[i]):
            if s != "" and d != "":
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

    return file_df, indices_df, rows, columns


def main(base_path, l=53, u=51, w=325, step=447, plot_show=False):

    file_df, indices_df, rows, columns = read_base_path(base_path)
    time_cos_sims = {}
    for i in range(len(file_df)):
        print(f"Processing {file_df.iloc[i].tif_path} at time {file_df.iloc[i].time}")
        encodings = crop_and_encode(
            file_df.iloc[i].tif_path,
            l=l,
            u=u,
            w=w,
            step=step,
            plot_show=plot_show,
            rows=rows,
            columns=columns,
        )
        time_cos_sim = compute_cos_sim(encodings, indices_df)
        time_cos_sims[file_df.iloc[i].time] = time_cos_sim
        time_cos_sim.to_csv(
            base_path / f"resnet50_cosine_similarity_{file_df.iloc[i].folder}.csv"
        )

    return time_cos_sims


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
    main(folder_path, plot_show=False)
