from mcw_utils.parsef import Parsef
import pandas as pd
from pathlib import Path
import argparse


def main(fname_sorted, fname_seq, ras, out_dir):
    if ras == "kras":
        isoform_list = [
            "ENST00000256078.8",
            "ENST00000311936.7",
            "ENST00000556131.1",
            "ENST00000557334.5",
        ]
    elif ras == "hras":
        isoform_list = [  #'ENST00000311189.7',
            #'ENST00000397594.5',
            #'ENST00000397596.6',
            #'ENST00000493230.5',
            "ENST00000451590.5",
            "ENST00000417302.5",
            "ENST00000462734.1",
            "ENST00000468682.2",
            "ENST00000478324.5",
            "ENST00000479482.1",
            #'ENST00000482021.1',
            "ENST00000369535.4",
        ]
        # isoform_list = [
        #    "ENST00000311189.7",
        #    "ENST00000397594.5",
        # ]
        # "ENST00000397596.6",
        # "ENST00000493230.5",
        # "ENST00000451590.5",
        # "ENST00000417302.5",
        #'ENST00000462734.1',
        #'ENST00000468682.2',
        #'ENST00000478324.5',
        #'ENST00000479482.1',
        #'ENST00000482021.1',
        #'ENST00000369535.4']
    else:
        raise ValueError("Unsupported RAS type. Use 'kras' or 'hras'.")
    isoforms = Parsef(
        bam_fname=fname_sorted, seq_fname=fname_seq, isoform_list=isoform_list
    )
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
    isoforms.measure_distributions()
    isoforms.fit_alignment_model()
    print("Simming experiment")
    isoforms.simulate_experiment(max_iter=5)
    if out_dir is not None:
        print("saving simulated data")
        isoforms.save_simulated_data(out_dir / f"{ras}_rsem_simulated_data.csv")
    print("Measuring isoform fractions")
    isoforms.measure_isoform_fractions(quantile=0.5)
    if out_dir is not None:
        print("Saving isoform fractions and effective lengths")
        isoforms.isoform_fractions.to_csv(
            out_dir / f"{ras}_rsem_isoform_fractions.csv", index=False
        )
        pd.DataFrame(
            isoforms.effective_lengths, index=["effective_length"]
        ).transpose().to_csv(out_dir / f"{ras}_rsem_effective_lengths.csv")


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
        "--ras",
        type=str,
        choices=["kras", "hras"],
        default="kras",
        help="Type of RAS isoform to measure.",
    )
    parser.add_argument(
        "--out_dir",
        type=str,
        default=None,
        help="Output directory to save results.",
    )

    args = parser.parse_args()
    main(args.fname_sorted, args.fname_seq, args.ras, args.out_dir)
