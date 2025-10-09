import pandas as pd
from pathlib import Path
from datetime import date


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
    CYBB = "ENSG00000165168"
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

    # RNA modifications - https://pmc.ncbi.nlm.nih.gov/articles/PMC10447785/pdf/43556_2023_Article_139.pdf fig1
    ## m1A
    ### writers
    TRMT10C = "ENSG00000174173"
    TRMT6 = "ENSG00000089195"
    TRMT61 = "ENSG00000166166"
    ### readers
    YTHD1 = "ENSG00000149658"
    YTHD3 = "ENSG00000185728"
    ### erasers
    ALKBH1 = "ENSG00000100601"
    ALKBH3 = "ENSG00000166199"
    ALKBH7 = "ENSG00000125652"

    ## m7G
    ### writers
    METTL1 = "ENSG00000037897"
    WBSCR22 = "ENSG00000071462"
    RNMT = "ENSG00000101654"
    ### readers
    ### erasers

    ## AtoI
    ### writers
    ADAR1 = "ENSG00000160710"
    ADAR2 = "ENSG00000197381"
    ADAT2 = "ENSG00000189007"
    ### readers
    ### erasers

    ## m6Am
    ### writers
    METTL4 = "ENSG00000101574"
    PCIF1 = "ENSG00000100982"
    ### readers
    ### erasers
    FTO = "ENSG00000140718"

    # psi
    ### writers
    PUS1 = "ENSG00000177192"
    PUS7L = "ENSG00000129317"
    RPUSD1 = "ENSG00000007376"
    DKC1 = "ENSG00000130826"
    ### readers
    ### erasers

    ## m6A
    ### writers
    METTL3 = "ENSG00000165819"
    METTL5 = "ENSG00000138382"
    METTL14 = "ENSG00000145388"
    METTL16 = "ENSG00000127804"
    ZCCHC4 = "ENSG00000168228"
    ### readers
    YTHDF1 = "ENSG00000149658"
    YTHDC1 = "ENSG00000275272"
    YTHDC2 = "ENSG00000047188"
    ### erasers
    ALKBH5 = "ENSG00000091542"

    # m5C
    ### writers
    NSUN1 = "ENSG00000111641"
    DNMT2 = "ENSG00000107614"
    ### readers
    # ALYREF = ""
    YTHDF2 = "ENSG00000198492"
    YBX1 = "ENSG00000065978"
    ### erasers
    TET1 = "ENSG00000138336"
    TET2 = "ENSG00000168769"

    # mcm5s2U
    ### writers
    ELP1 = "ENSG00000070061"
    ELP3 = "ENSG00000134014"
    CTU1 = "ENSG00000142544"
    CTU2 = "ENSG00000174177"
    ### readers
    ### erasers

    # mcm5U
    ### writers
    ALKBH8 = "ENSG00000137760"
    ### readers
    ### erasers

    # ac4C
    ### writers
    NAT10 = "ENSG00000135372"
    ### readers
    ### erasers
    SIRT7 = "ENSG00000187531"

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
        "CYBB": CYBB,  # CYBB is an alias for NOX2
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

    m1a = {
        "TRMT10C": TRMT10C,
        "TRMT6": TRMT6,
        "TRMT61": TRMT61,
        "YTHD1": YTHD1,
        "YTHD3": YTHD3,
        "ALKBH1": ALKBH1,
        "ALKBH3": ALKBH3,
        "ALKBH7": ALKBH7,
    }

    m7g = {
        "METTL1": METTL1,
        "WBSCR22": WBSCR22,
        "RNMT": RNMT,
    }

    AtoI = {
        "ADAR1": ADAR1,
        "ADAR2": ADAR2,
        "ADAT2": ADAT2,
    }

    m6Am = {
        "METTL4": METTL4,
        "PCIF1": PCIF1,
        "FTO": FTO,
    }

    psi = {
        "PUS1": PUS1,
        "PUS7L": PUS7L,
        "RPUSD1": RPUSD1,
        "DKC1": DKC1,
    }

    m6a = {
        "METTL3": METTL3,
        "METTL5": METTL5,
        "METTL14": METTL14,
        "METTL16": METTL16,
        "ZCCHC4": ZCCHC4,
        "YTHDF1": YTHDF1,
        "YTHDC1": YTHDC1,
        "YTHDC2": YTHDC2,
        "ALKBH5": ALKBH5,
    }

    m5c = {
        "NSUN1": NSUN1,
        "DNMT2": DNMT2,
        "YTHDF2": YTHDF2,
        "YBX1": YBX1,
        "TET1": TET1,
        "TET2": TET2,
    }

    mcm5s2u = {
        "ELP1": ELP1,
        "ELP3": ELP3,
        "CTU1": CTU1,
        "CTU2": CTU2,
    }

    mcm5u = {
        "ALKBH8": ALKBH8,
    }

    ac4c = {
        "NAT10": NAT10,
        "SIRT7": SIRT7,
    }

    data = {}
    data["accession_id"] = []
    # data = init_data(up, data)
    # data = init_data(down, data)
    # data = init_data(housekeeping, data)
    # data = init_data(sweeny, data)
    # data = init_data(adriano, data)
    # data = init_data(tommy, data)
    data = init_data(m1a, data)
    data = init_data(m7g, data)
    data = init_data(AtoI, data)
    # data = init_data(m6Am, data)
    data = init_data(psi, data)
    # data = init_data(m6a, data)
    data = init_data(m5c, data)
    data = init_data(mcm5s2u, data)
    data = init_data(mcm5u, data)
    data = init_data(ac4c, data)

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
                            # data = pull_tpm(up, dfi, data)
                            # data = pull_tpm(down, dfi, data)
                            # data = pull_tpm(housekeeping, dfi, data)
                            # data = pull_tpm(sweeny, dfi, data)
                            # data = pull_tpm(adriano, dfi, data)
                            # data = pull_tpm(tommy, dfi, data)
                            data = pull_tpm(m1a, dfi, data)
                            data = pull_tpm(m7g, dfi, data)
                            data = pull_tpm(AtoI, dfi, data)
                            # data = pull_tpm(m6Am, dfi, data)
                            data = pull_tpm(psi, dfi, data)
                            # data = pull_tpm(m6a, dfi, data)
                            data = pull_tpm(m5c, dfi, data)
                            data = pull_tpm(mcm5s2u, dfi, data)
                            data = pull_tpm(mcm5u, dfi, data)
                            data = pull_tpm(ac4c, dfi, data)

    df = pd.DataFrame(data)
    fname = f"gene_expression_results_{date.today()}.csv"
    df.to_csv(fname, index=False)
    print(f"Gene expression results saved to {fname}")


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
