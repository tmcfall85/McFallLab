from mcw_utils.parsef import Parsef
import pandas as pd
from pathlib import Path
import argparse


def main(fname_sorted, fname_seq, out_dir):
    isoform_list = ["ENST00000256078.8"]
    isoforms = Parsef(
        bam_fname=fname_sorted, seq_fname=fname_seq, isoform_list=isoform_list
    )

    print(f"Creating output directory: {out_dir}")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print("Measuring distributions")
    isoforms.measure_distributions()
    left_hist_dist = isoforms.left_start_distributions[isoform_list[0]]
    df = pd.DataFrame(
        {"count": list(left_hist_dist[0]) + [0], "edge": left_hist_dist[1]}
    ).to_csv(out_dir / f"nras_left_start_distribution.csv", index=False)

    right_hist_dist = isoforms.right_start_distributions[isoform_list[0]]
    df = pd.DataFrame(
        {"count": list(right_hist_dist[0]) + [0], "edge": right_hist_dist[1]}
    ).to_csv(out_dir / f"nras_right_start_distribution.csv", index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Measure RAS isoform stoichiometry.")
    parser.add_argument(
        "--fname_sorted",
        type=str,
        help="Path to the sorted BAM file.",
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
    main(args.fname_sorted, args.fname_seq, args.out_dir)
