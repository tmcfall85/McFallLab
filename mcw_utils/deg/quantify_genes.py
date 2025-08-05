import pandas as pd
from pathlib import Path
from datetime import date


def main(folder_path):

    data = {}
    data["accession_id"] = []
    needs_init = True
    genes = []
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
                            data["accession_id"].append(acc_num)

                            if needs_init:
                                genes = list(set(dfi.gene_id))
                                for gene in genes:
                                    data[gene] = []
                                needs_init = False

                            for gene in genes:
                                gene_results = dfi[dfi.gene_id.str.startswith(gene)]
                                if len(gene_results) > 0:
                                    data[gene].append(gene_results.iloc[0].TPM)
                                else:
                                    data[gene].append(0)

    df = pd.DataFrame(data)
    out_file = f"deg_gene_expression_results_{date.today()}.csv"
    df.to_csv(out_file, index=False)
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
