from mcw_utils.parsef import Parsef
import pandas as pd
from pathlib import Path
import argparse

target_directory = Path("path/to/your/directory")  # Replace with your directory

for item in target_directory.iterdir():
    if item.is_dir():
        print(f"Found folder: {item}")


def main(target_dir, fname_seq, out_dir):
    isoform_list = ["ENST00000256078.8"]
    for item in target_directory.iterdir():
        if item.is_dir():
            print(f"Found folder: {item}")
            fname_sorted = (
                item / "star_output" / "Aligned.toTranscriptome.sorted.out.bam"
            )
            isoforms = Parsef(
                bam_fname=fname_sorted, seq_fname=fname_seq, isoform_list=isoform_list
            )

            print(f"Creating output directory: {out_dir}")
            out_dir = item / out_dir
            out_dir.mkdir(parents=True, exist_ok=True)
            print("Measuring distributions")
            isoforms.measure_distributions()
            left_hist_dist = isoforms.left_start_distributions[isoform_list[0]]
            df = pd.DataFrame(
                {"count": list(left_hist_dist[0]) + [0], "edge": left_hist_dist[1]}
            ).to_csv(out_dir / f"kras_left_start_distribution.csv", index=False)

            right_hist_dist = isoforms.right_start_distributions[isoform_list[0]]
            df = pd.DataFrame(
                {"count": list(right_hist_dist[0]) + [0], "edge": right_hist_dist[1]}
            ).to_csv(out_dir / f"kras_right_start_distribution.csv", index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Measure RAS isoform stoichiometry.")
    parser.add_argument(
        "--target_dir",
        type=str,
        help="Path to the directory containing sorted BAM files.",
    )
    parser.add_argument(
        "--fname_seq",
        type=str,
        help="Path to the sequence file.",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        help="Output directory to save results.",
    )

    args = parser.parse_args()
    main(args.target_dir, args.fname_seq, args.out_dir)
