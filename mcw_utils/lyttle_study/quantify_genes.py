import pandas as pd
from pathlib import Path


def main(folder_path):
    # up
    TPO = "ENSG00000115705"
    MMP1 = "ENSG00000196611"
    CCL17 = "ENSG00000102970"
    CCL23 = "ENSG00000274736"
    DKK1 = "ENSG00000107984"
    SERPINE1 = "ENSG00000106366"
    TNFSF12 = "ENSG00000239697"
    SOST = "ENSG00000167941"
    ICAM1 = "ENSG00000090339"
    SPARC = "ENSG00000113140"
    TIMP1 = "ENSG00000102265"

    # down
    BGLAP = "ENSG00000242252"
    IL16 = "ENSG00000172349"
    IL1A = "ENSG00000115008"
    IL9 = "ENSG00000145839"
    CCL3 = "ENSG00000278567"
    IL35 = "ENSG00000138378"
    IL23R = "ENSG00000162594"
    IL23A = "ENSG00000110944"
    IL17A = "ENSG00000112115"
    IL1B = "ENSG00000125538"

    up = {
        "TPO": TPO,
        "MMP1": MMP1,
        "CCL17": CCL17,
        "CCL23": CCL23,
        "DKK1": DKK1,
        "SERPINE1": SERPINE1,
        "TNFSF12": TNFSF12,
        "SOST": SOST,
        "ICAM1": ICAM1,
        "SPARC": SPARC,
        "TIMP1": TIMP1,
    }

    down = {
        "BGLAP": BGLAP,
        "IL16": IL16,
        "IL1A": IL1A,
        "IL9": IL9,
        "CCL3": CCL3,
        "IL35": IL35,
        "IL23R": IL23R,
        "IL23A": IL23A,
        "IL1LA": IL17A,
        "IL1B": IL1B,
    }

    data = {}
    data["accession_id"] = []
    for gene in up.keys():
        col_name = f"{gene}_TPM"
        data[col_name] = []
    for gene in down.keys():
        col_name = f"{gene}_TPM"
        data[col_name] = []
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
                            print(dfi.iloc[0])

                            for gene in up.keys():
                                gene_results = dfi[dfi.gene_id.str.startswith(up[gene])]
                                col_name = f"{gene}_TPM"
                                if len(gene_results) > 0:

                                    data[col_name].append(gene_results.iloc[0].TPM)
                                else:
                                    data[col_name].append(0)

                            for gene in down.keys():
                                gene_results = dfi[
                                    dfi.gene_id.str.startswith(down[gene])
                                ]
                                col_name = f"{gene}_TPM"
                                if len(gene_results) > 0:

                                    data[col_name].append(gene_results.iloc[0].TPM)
                                else:
                                    data[col_name].append(0)

    df = pd.DataFrame(data)
    df.to_csv("gene_expression_results.csv", index=False)
    print("Gene expression results saved to gene_expression_results.csv")


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
