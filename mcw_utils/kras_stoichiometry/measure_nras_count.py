from mcw_utils.parsef import Parsef
import pandas as pd
from pathlib import Path
import argparse


def main(target_dir, fname_seq, out_dir):
    target_dir = Path(target_dir)
    isoform_list = ["ENST00000369535.4"]
    left_count = {}
    right_count = {}
    left_edge = {}
    right_edge = {}
    for item in target_dir.iterdir():
        if item.is_dir():
            acc_num = item.stem
            print(f"Found folder: {item}")
            fname_sorted = (
                item / "star_output" / "Aligned.toTranscriptome.sorted.out.bam"
            )
            try:
                isoforms = Parsef(
                    bam_fname=fname_sorted,
                    seq_fname=fname_seq,
                    isoform_list=isoform_list,
                )

                print("Measuring distributions")
                isoforms.measure_distributions()
                left_hist_dist = isoforms.left_start_distributions[isoform_list[0]]
                right_hist_dist = isoforms.right_start_distributions[isoform_list[0]]
                left_count[acc_num] = list(left_hist_dist[0])
                right_count[acc_num] = list(right_hist_dist[0])
                left_edge[acc_num] = list(left_hist_dist[1])
                right_edge[acc_num] = list(right_hist_dist[1])
            except Exception as e:
                print(f"Error processing {acc_num}: {e}")
                continue

    print(f"Creating output directory: {out_dir}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(left_count).to_csv(
        out_dir / f"nras_left_start_distribution.csv", index=False
    )

    pd.DataFrame(right_count).to_csv(
        out_dir / f"nras_right_start_distribution.csv", index=False
    )
    pd.DataFrame(left_edge).to_csv(out_dir / f"nras_left_start_edge.csv", index=False)
    pd.DataFrame(right_edge).to_csv(out_dir / f"nras_right_start_edge.csv", index=False)


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
