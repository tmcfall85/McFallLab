import anndata
import pandas as pd


def remove_ver(row):
    return row["Unnamed: 0"].split(".")[0]


def is_panc(row):
    return row.specimen_sample_site.startswith("Panc")


deg_data = pd.read_csv(
    "/home/msochor/repos/McFallLab/mcw_utils/deg/deg_results_with_metadata_2025-08-07.csv"
)
adata_df = anndata.read_h5ad(
    "/mnt/c/Users/msochor/Downloads/dominguez_conde_immune_tissue_two_donors.h5ad"
).var
big_data = pd.read_csv(
    "/home/msochor/repos/McFallLab/mcw_utils/deg/deg_gene_expression_results_merged_to_outcomes.csv",
    index_col="accession_id",
).drop(columns="Unnamed: 0")
linker_df = pd.read_csv(
    "/home/msochor/repos/scratch/output_8_4_2025/tempus_json_pdac_rcc_merged_dna_tall.csv"
)
unfiltered_df = pd.read_stata(
    "/mnt/c/Users/msochor/Downloads/PDAC_SOTB_2025_03_04_unfiltered.dta"
)

deg_data.sort_values(by="padj", ascending=True, inplace=True)
top_1000_sig_genes = deg_data.head(1000).copy()
top_1000_sig_genes["gene_code"] = top_1000_sig_genes.apply(remove_ver, axis=1)
top_500_sig_genes = deg_data.head(500).copy()
top_500_sig_genes["gene_code"] = top_500_sig_genes.apply(remove_ver, axis=1)


linker_df.specimen_sample_site = linker_df.specimen_sample_site.fillna("unknown")
linker_df["is_panc"] = linker_df.apply(is_panc, axis=1)
linker_rna_df = linker_df[linker_df["Has RNA Folder"] == "Yes"]
unique_panc_linker_rna_df = (
    linker_rna_df.groupby(["mrn", "accession_id"]).first().reset_index()
)

unique_unfiltered = (
    unfiltered_df[["mrn", "pfsmofromdx"]].groupby("mrn").first().reset_index()
)

big_data_w_mo = big_data.merge(
    unique_panc_linker_rna_df[
        ["accession_id", "mrn", "specimen_sample_site", "is_panc"]
    ],
    on="accession_id",
    how="left",
).merge(unique_unfiltered, on="mrn", how="left")

big_data_w_mo.specimen_sample_site = big_data_w_mo.specimen_sample_site.fillna(
    "unknown"
)
big_data_w_mo_T = big_data_w_mo.set_index("accession_id").T

big_data_w_mo.to_csv(
    "/mnt/c/Users/msochor/Downloads/big_data_with_labels.csv", index=False
)

gene_names = []
top_1000_sig = []
top_500_sig = []
for i in big_data_w_mo_T.index:
    if i.find("ENSG") > -1:
        gene_code = i.split(".")[0]
        if gene_code in top_1000_sig_genes["gene_code"].values:
            top_1000_sig.append(True)
        else:
            top_1000_sig.append(False)
        if gene_code in top_500_sig_genes["gene_code"].values:
            top_500_sig.append(True)
        else:
            top_500_sig.append(False)
        gene_name = adata_df[adata_df.ensembl_id == gene_code].gene_name
        if len(gene_name) > 0:
            gene_names.append(gene_name.values[0].split("_ENS")[0])
        else:
            gene_names.append(None)
        # else:
        #    gene_names.append(None)
    else:
        gene_names.append(None)
        top_1000_sig.append(False)
        top_500_sig.append(False)

big_data_w_mo_T["gene_name"] = gene_names
big_data_w_mo_T["top_1000_sig"] = top_1000_sig
big_data_w_mo_T["top_500_sig"] = top_500_sig
big_data_genes = big_data_w_mo_T[big_data_w_mo_T.gene_name.notnull()].copy()
big_data_genes.set_index("gene_name", inplace=True)

big_data_genes.to_csv(
    "/mnt/c/Users/msochor/Downloads/big_data_for_c2s_modeling.csv", index=True
)
