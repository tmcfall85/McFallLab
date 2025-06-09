from mcw_utils.parsef import Parsef
import pandas as pd
from pathlib import Path

fname_sorted = (
    "/mnt/c/Users/msochor/Downloads/tlsort3/TL-20-8D3FFE_rsem.transcript.sorted.bam"
)
fname_seq = "/mnt/c/Users/msochor/Downloads/rsem.seq"


def main(fname_sorted=fname_sorted, fname_seq=fname_seq, ras="kras", out_dir=None):
    if ras == "kras":
        isoform_list = [
            "ENST00000256078.8",
            "ENST00000311936.7",
            "ENST00000556131.1",
            "ENST00000557334.5",
        ]
    elif ras == "hras":
        isoform_list = [
            "ENST00000311189.7",
            "ENST00000397594.5",
        ]
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
    isoforms.measure_distributions()
    isoforms.fit_alignment_model()
    isoforms.simulate_experiment(max_iter=5)
    isoforms.measure_isoform_fractions(quantile=0.5)
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        isoforms.save_simulated_data(out_dir / f"{ras}_rsem_simulated_data.csv")
        isoforms.isoform_fractions.to_csv(
            out_dir / f"{ras}_rsem_isoform_fractions.csv", index=False
        )
        pd.DataFrame(isoforms.effective_lengths).to_csv(
            out_dir / f"{ras}_rsem_effective_lengths.csv", index=False
        )
