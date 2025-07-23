import pandas as pd
from pathlib import Path


def init_data(gene_dict, data):
    for gene in gene_dict.keys():
        col_name = f"{gene}_TPM"
        data[col_name] = []
    return data


def pull_tpm(gene_dict, dfi, data):
    for gene in gene_dict.keys():
        gene_results = dfi[dfi.gene_id.str.startswith(gene_dict[gene])]
        col_name = f"{gene}_TPM"
        if len(gene_results) > 0:

            data[col_name].append(gene_results.iloc[0].TPM)
        else:
            data[col_name].append(0)
    return data


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

    # housekeeping
    ACTB = "ENSG00000075624"
    GAPDH = "ENSG00000111640"
    RPLP0 = "ENSG00000089157"
    PGK1 = "ENSG00000102144"
    PPIA = "ENSG00000196262"
    RPL13A = "ENSG00000142541"
    RLA0 = "ENSG00000089157"

    # Sweeny lab
    NOX5 = "ENSG00000255346"
    NOX4 = "ENSG00000086991"
    NOX3 = "ENSG00000074771"
    NOX1 = "ENSG00000007952"

    # Adriano Lab
    ARRB1 = "ENSG00000137486"
    ARRB2 = "ENSG00000141480"
    ARBK1 = "ENSG00000173020"
    GRK4 = "ENSG00000125388"
    GRK1 = "ENSG00000288263"
    GRK5 = "ENSG00000198873"
    ARBK2 = "ENSG00000100077"
    GRK6 = "ENSG00000198055"
    GRK7 = "ENSG00000114124"
    ADRB1 = "ENSG00000043591"
    ADRB2 = "ENSG00000169252"
    ADRB3 = "ENSG00000188778"

    # Tommy Lab
    KRAS = "ENSG00000133703.11"
    HRAS = "ENSG00000174775.16"
    NRAS = "ENSG00000213281.4"

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

    housekeeping = {
        "ACTB": ACTB,
        "GAPDH": GAPDH,
        "RPLP0": RPLP0,
        "PGK1": PGK1,
        "PPIA": PPIA,
        "RPL13A": RPL13A,
        "RLA0": RLA0,
    }

    sweeny = {
        "NOX5": NOX5,
        "NOX4": NOX4,
        "NOX3": NOX3,
        "NOX1": NOX1,
    }

    adriano = {
        "ARRB1": ARRB1,
        "ARRB2": ARRB2,
        "ARBK1": ARBK1,
        "GRK4": GRK4,
        "GRK1": GRK1,
        "GRK5": GRK5,
        "ARBK2": ARBK2,
        "GRK6": GRK6,
        "GRK7": GRK7,
        "ADRB1": ADRB1,
        "ADRB2": ADRB2,
        "ADRB3": ADRB3,
    }

    tommy = {
        "KRAS": KRAS,
        "HRAS": HRAS,
        "NRAS": NRAS,
    }

    data = {}
    data["accession_id"] = []
    data = init_data(up, data)
    data = init_data(down, data)
    data = init_data(housekeeping, data)
    data = init_data(sweeny, data)
    data = init_data(adriano, data)
    data = init_data(tommy, data)

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
                            data = pull_tpm(up, dfi, data)
                            data = pull_tpm(down, dfi, data)
                            data = pull_tpm(housekeeping, dfi, data)
                            data = pull_tpm(sweeny, dfi, data)
                            data = pull_tpm(adriano, dfi, data)
                            data = pull_tpm(tommy, dfi, data)

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
