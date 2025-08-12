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


def measure_vec(st):
    """
    Measure the vector representations of images across time points.

    Note: different images and contrasts make this approach impractical and maybe not correct.
    I keep this function here for reference, if at some later point I want to see this approach.
    """
    encodings = []
    fig, axs = plt.subplots(1, 6)
    for i in range(6):
        img_file = f"/mnt/c/Users/msochor/Downloads/6.19.25_before drug/0005587_01/{st}{i}.png"  # Replace with your image URL or path
        image = Image.open(img_file).convert("RGB")
        axs[i].imshow(image)
        image_processor = AutoImageProcessor.from_pretrained("microsoft/resnet-50")
        model = ResNetModel.from_pretrained("microsoft/resnet-50")

        encoding = []
        for i in range(4):
            image = image.transpose(Image.ROTATE_90)
            inputs = image_processor(image, return_tensors="pt")
            with torch.no_grad():
                outputs = model(**inputs)
                encoding.append(outputs.pooler_output.squeeze())

        encodings.append(encoding)
    stacked_tensors = torch.stack(encodings[0])
    average_zero_tensor = torch.mean(stacked_tensors, dim=0)
    # return encodings
    plt.show()
    mean_cos_sims = []
    std_cos_sims = []
    for i in range(6):
        cos_sims = []
        for encoding in encodings[i]:
            # for rotation in encoding:
            cos_sims.append(cos_sim(average_zero_tensor, encoding))
        mean_cos_sims.append(np.mean(cos_sims))
        std_cos_sims.append(np.std(cos_sims))

    print(mean_cos_sims)
    print(std_cos_sims)
    return mean_cos_sims, std_cos_sims


def measure_vec_timepoint(folder_path, tp, replicates, plot_show=False):
    """
    Measure the vector representations of images at each time point relative to vehicle.

    This approach is superior because the vehicle degrades with time as well,
    so in effect this measures how much additional degradation occurs with drug treatment.
    """
    encodings = []
    fig, axs = plt.subplots(3, 2)
    prefixes = ["vv", "rr", "tt"]
    image_processor = AutoImageProcessor.from_pretrained("microsoft/resnet-50")
    model = ResNetModel.from_pretrained("microsoft/resnet-50")
    for i, (drug, replicate) in enumerate(zip(prefixes, replicates)):
        encoding = []
        for j in range(2):
            if j + 1 in replicate:
                img_file = (
                    folder_path / f"{drug}{j+1}{tp}.png"
                )  # Replace with your image URL or path
                image = Image.open(img_file).convert("RGB")

                axs[i, j].imshow(image)

                for k in range(4):
                    image = image.transpose(Image.ROTATE_90)
                    inputs = image_processor(image, return_tensors="pt")
                    with torch.no_grad():
                        outputs = model(**inputs)
                        encoding.append(outputs.pooler_output.squeeze())

        encodings.append(encoding)
    stacked_tensors = torch.stack(encodings[0])
    average_zero_tensor = torch.mean(stacked_tensors, dim=0)
    if plot_show:
        plt.show()
    mean_cos_sims = []
    std_cos_sims = []
    for i in range(3):
        cos_sims = []
        for encoding in encodings[i]:
            # for rotation in encoding:
            cos_sims.append(cos_sim(average_zero_tensor, encoding))
        mean_cos_sims.append(np.mean(cos_sims))
        std_cos_sims.append(np.std(cos_sims))
    print(mean_cos_sims)
    print(std_cos_sims)
    return mean_cos_sims, std_cos_sims


def main(folder_path, replicates, plot_show=False):
    categories = ["vehicle", "RMC", "TRM"]
    t0_mean, t0_std = measure_vec_timepoint(
        folder_path, "0", replicates=replicates, plot_show=plot_show
    )
    t1_mean, t1_std = measure_vec_timepoint(
        folder_path, "1", replicates=replicates, plot_show=plot_show
    )
    t2_mean, t2_std = measure_vec_timepoint(
        folder_path, "2", replicates=replicates, plot_show=plot_show
    )
    t3_mean, t3_std = measure_vec_timepoint(
        folder_path, "3", replicates=replicates, plot_show=plot_show
    )
    t4_mean, t4_std = measure_vec_timepoint(
        folder_path, "4", replicates=replicates, plot_show=plot_show
    )
    t5_mean, t5_std = measure_vec_timepoint(
        folder_path, "5", replicates=replicates, plot_show=plot_show
    )
    df_mean = pd.DataFrame(
        {
            "0": t0_mean,
            "0.5": t1_mean,
            "12": t2_mean,
            "48": t3_mean,
            "72": t4_mean,
            "96": t5_mean,
        },
        index=categories,
    ).transpose()
    df_std = pd.DataFrame(
        {
            "0": t0_std,
            "0.5": t1_std,
            "12": t2_std,
            "48": t3_std,
            "72": t4_std,
            "96": t5_std,
        },
        index=categories,
    ).transpose()
    df_mean["vehicle_std"] = df_std["vehicle"]
    df_mean["RMC_std"] = df_std["RMC"]
    df_mean["TRM_std"] = df_std["TRM"]
    s = []
    for r, c in zip(replicates, categories):
        rs = [str(x) for x in r]
        s.append(c + "".join(rs))
    fname_out = "_".join(s) + ".csv"
    df_mean.to_csv(fname_out)
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
        "replicates",
        type=str,
        help="List of lists containing the replicate numbers.",
        default="[[1,2],[1,2],[1,2]]",
    )
    args = parser.parse_args()
    folder_path = Path(args.folder_path)
    replicates = eval(args.replicates)
    main(folder_path, replicates)
