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


def compute_cos_sim(average_zero_tensor, encodings0, encodings1):
    stack = []
    for t in encodings0[0]:
        stack.append(t)

    stacked_tensors = torch.stack(stack)
    average_zero_tensor = torch.mean(stacked_tensors, dim=0)
    mean_cos_sims = []
    std_cos_sims = []
    for i in range(len(encodings0)):
        cos_sims = []
        for encoding in encodings0[i]:
            # for rotation in encoding:
            cos_sims.append(cos_sim(average_zero_tensor, encoding))
        mean_cos_sims.append(np.mean(cos_sims))
        std_cos_sims.append(np.std(cos_sims))
    for i in range(len(encodings1)):
        cos_sims = []
        for encoding in encodings1[i]:
            # for rotation in encoding:
            cos_sims.append(cos_sim(average_zero_tensor, encoding))
        mean_cos_sims.append(np.mean(cos_sims))
        std_cos_sims.append(np.std(cos_sims))

    print(mean_cos_sims)
    print(std_cos_sims)
    return mean_cos_sims, std_cos_sims


def measure_vec_timepoint(img_file, l=560, u=560, w=200, step=448, plot_show=False):
    """
    Measure the vector representations of images at each time point relative to vehicle.

    This is for the experiment data Raven emailed to be on 7-16-25 which was a 'previous PDAC spheroid experiment'.
    These images were cropped from a larger image, but the code here still autocrops wells from the
    larger image.
    """
    encodings = []
    fig, axs = plt.subplots(6, 10)

    # img_file = folder_path / f"{st}.png"  # Replace with your image URL or path
    image = Image.open(img_file).convert("RGB")
    image_processor = AutoImageProcessor.from_pretrained("microsoft/resnet-50")
    model = ResNetModel.from_pretrained("microsoft/resnet-50")
    for j in range(10):
        encoding = []
        for i in range(6):
            cropped = image.crop(
                (l + j * step, u + i * step, l + j * step + w, u + i * step + w)
            )
            if plot_show:
                axs[i, j].imshow(cropped)

            for i in range(4):
                cropped = cropped.transpose(Image.ROTATE_90)
                inputs = image_processor(cropped, return_tensors="pt")
                with torch.no_grad():
                    outputs = model(**inputs)
                    encoding.append(outputs.pooler_output.squeeze())

        encodings.append(encoding)
    return encodings


def main(folder_path, fname_out, plot_show=False):

    t00 = measure_vec_timepoint(folder_path, "00", plot_show=plot_show)
    t01 = measure_vec_timepoint(folder_path, "01", plot_show=plot_show)
    t10 = measure_vec_timepoint(folder_path, "10", plot_show=plot_show)
    t11 = measure_vec_timepoint(folder_path, "11", plot_show=plot_show)
    t20 = measure_vec_timepoint(folder_path, "20", step=112, plot_show=plot_show)
    t21 = measure_vec_timepoint(folder_path, "21", step=112, plot_show=plot_show)
    t30 = measure_vec_timepoint(folder_path, "30", step=112, plot_show=plot_show)
    t31 = measure_vec_timepoint(folder_path, "31", step=112, plot_show=plot_show)
    stack = []
    for t in t00:
        for tt in t:
            stack.append(tt)
    for t in t01:
        for tt in t:
            stack.append(tt)

    stacked_tensors = torch.stack(stack)
    average_zero_tensor = torch.mean(stacked_tensors, dim=0)

    time0_mean, time0_std = compute_cos_sim(average_zero_tensor, t00, t01)
    time1_mean, time1_std = compute_cos_sim(average_zero_tensor, t10, t11)
    time2_mean, time2_std = compute_cos_sim(average_zero_tensor, t20, t21)
    time3_mean, time3_std = compute_cos_sim(average_zero_tensor, t31, t31)
    df_mean = pd.DataFrame(
        {
            "before": time0_mean,
            "after": time1_mean,
            "48": time2_mean,
            "72": time3_mean,
        },
        index=[
            "vehicle",
            "0.1nM",
            "0.5nM",
            "2nM",
            "5nM",
            "10nM",
            "15nM",
            "20nM",
            "25nM",
            "30nM",
        ],
    ).transpose()
    df_std = pd.DataFrame(
        {
            "before": time0_std,
            "after": time1_std,
            "48": time2_std,
            "72": time3_std,
        },
        index=[
            "vehicle",
            "0.1nM",
            "0.5nM",
            "2nM",
            "5nM",
            "10nM",
            "15nM",
            "20nM",
            "25nM",
            "30nM",
        ],
    ).transpose()
    df_mean["vehicle_std"] = df_std["vehicle"]
    df_mean["0.1nM_std"] = df_std["0.1nM"]
    df_mean["0.5nM_std"] = df_std["0.5nM"]
    df_mean["2nM_std"] = df_std["2nM"]
    df_mean["5nM_std"] = df_std["5nM"]
    df_mean["10nM_std"] = df_std["10nM"]
    df_mean["15nM_std"] = df_std["15nM"]
    df_mean["20nM_std"] = df_std["20nM"]
    df_mean["25nM_std"] = df_std["25nM"]
    df_mean["30nM_std"] = df_std["30nM"]
    df_mean.transpose().to_csv(fname_out)
    print(f"Results saved to {fname_out}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Measure donuts in file path.")
    parser.add_argument(
        "folder_path",
        type=str,
        help="Path to the folder containing cropped donut images.",
    )
    parser.add_argument(
        "fname_out",
        type=str,
        help="Filename to save results to.",
    )
    args = parser.parse_args()
    folder_path = Path(args.folder_path)
    main(folder_path, args.fname_out, plot_show=False)
