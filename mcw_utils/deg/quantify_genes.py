import pandas as pd
from pathlib import Path
from datetime import date


def main(folder_path):

    data = []
    needs_init = True
    search_dir = Path(folder_path)
    rsem_name = "Aligned.toTranscriptome.out_rsem.genes.results"
    for acc_dir in search_dir.iterdir():
        if acc_dir.is_dir():
            acc_num = acc_dir.stem
            for sub_folder in acc_dir.iterdir():
                if sub_folder.is_dir() and sub_folder.name == "rsem_output":
                    for fname in sub_folder.iterdir():
                        if fname.is_file() and fname.name == rsem_name:
                            print("File to be processed:", fname)

                            dfi = pd.read_csv(fname, sep="\t")
                            dfi = dfi[dfi.gene_id.str.startswith("ENSG0")]
                            dfi_t = dfi.set_index("gene_id")[["FPKM", "TPM"]].T
                            dfi_t["accession_id"] = acc_num
                            data.append(dfi_t)

    df = pd.concat(data).reset_index()
    df = df[df["index"] == "TPM"]
    df.drop(columns=["index"], inplace=True)
    df.fillna(0, inplace=True)
    df.set_index("accession_id", inplace=True)
    out_file = f"deg_gene_expression_results_{date.today()}.csv"
    df.to_csv(out_file)
    print(f"Gene expression results saved to {out_file}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract TPM values for specific genes from RSEM output."
    )
    parser.add_argument(
        "folder_path",
        type=str,
        help="Path to the folder containing subfolders with RSEM output.",
    )

    args = parser.parse_args()
    main(args.folder_path)
