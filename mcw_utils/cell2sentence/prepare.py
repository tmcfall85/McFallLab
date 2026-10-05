import anndata
import pandas as pd
from datetime import datetime
from pathlib import Path
import json
import sys


def remove_ver(row):
    return row["gene_id"].split(".")[0]


def main(deg_folder, top_n_sig):

    metadata = {}
    for file_name in deg_folder.iterdir():
        if file_name.is_file():
            if file_name.name == "metadata.json":
                with open(file_name, "r") as fp:
                    metadata = json.load(fp)
                print("Deg metadata:")
                print(metadata)
    results_dir = Path(f"prepare/{str(datetime.now()).replace(' ','_')}")
    results_dir.mkdir(exist_ok=True)
    for file_path in deg_folder.iterdir():
        if file_path.is_file():
            if file_path.name.startswith("deg_results"):
                print(f"Preparing deg file: {file_path.name}")
                prepare_metadata = prepare(file_path, top_n_sig, results_dir)

                metadata["deg_results_file"] = str(file_path)
                metadata["top_n_sig"] = top_n_sig
                for key in prepare_metadata.keys():
                    metadata[key] = prepare_metadata[key]
                with open(results_dir / f"{file_path.stem}_metadata.json", "w") as fp:
                    json.dump(metadata, fp)


def prepare(deg_file, top_n_sig, results_dir):
    deg_data = pd.read_csv(deg_file)
    deg_data.rename({"Unnamed: 0": "gene_id"}, axis=1, inplace=True)
    adata_file = (
        "/mnt/c/Users/msochor/Downloads/dominguez_conde_immune_tissue_two_donors.h5ad"
    )
    if Path(adata_file).is_file():
        adata_df = anndata.read_h5ad(adata_file).var
    else:
        adata_df = pd.DataFrame({'ensembl_id': [], 'gene_name': []})
    train_data = pd.read_csv(deg_file.parent / "train.csv")
    train_data_T = train_data.set_index("accession_id").T

    deg_data.sort_values(by="padj", ascending=True, inplace=True)

    holdout_file = deg_file.parent / "holdout.csv"
    metadata = {}
    train_annotated_filenames = []
    holdout_annotated_filenames = []
    for n in top_n_sig:
        top_sig_genes = deg_data.head(int(n)).copy()
        top_sig_genes["gene_code"] = top_sig_genes.apply(remove_ver, axis=1)

        train_annotated_T = annotate(train_data_T, adata_df, top_sig_genes)

        train_annotated = train_annotated_T.T
        train_merge_cols = [
            "accession_id",
            "pfsmofromdx",
            "osmofromdx",
            "recurrence_time_sur",
        ]
        train_annotated_with_outcomes = train_annotated.merge(
            train_data[train_merge_cols], on="accession_id", how="left"
        )
        train_annotated_filename = (
            results_dir / f"{deg_file.stem}_top_{n}_sig_train_for_c2s_modeling.csv"
        )
        train_annotated_with_outcomes.to_csv(
            train_annotated_filename,
            index=False,
        )
        train_annotated_filenames.append(str(train_annotated_filename))

        if holdout_file.is_file():
            holdout_data = pd.read_csv(holdout_file)
            holdout_data_T = holdout_data.set_index("accession_id").T
            holdout_annotated_T = annotate(holdout_data_T, adata_df, top_sig_genes)
            holdout_annotated = holdout_annotated_T.T
            holdout_merge_cols = [
                "accession_id",
                "pfsmofromdx",
                "osmofromdx",
                "recurrence_time_sur",
            ]
            holdout_annotated_with_outcomes = holdout_annotated.merge(
                holdout_data[holdout_merge_cols], on="accession_id", how="left"
            )
            holdout_annotated_filename = (
                results_dir
                / f"{deg_file.stem}_top_{n}_sig_holdout_for_c2s_modeling.csv"
            )
            holdout_annotated_with_outcomes.to_csv(
                holdout_annotated_filename,
                index=False,
            )
            holdout_annotated_filenames.append(str(holdout_annotated_filename))
    metadata["train_annotated_filenames"] = train_annotated_filenames
    metadata["holdout_annotated_filenames"] = holdout_annotated_filenames
    return metadata


def annotate(data_T, adata_df, top_sig_genes):
    gene_names = []
    top_sig = []
    for i in data_T.index:
        if i.find("ENSG") > -1:
            gene_code = i.split(".")[0]
            if gene_code in top_sig_genes["gene_code"].values:
                top_sig.append(True)
            else:
                top_sig.append(False)
            gene_name = adata_df[adata_df.ensembl_id == gene_code].gene_name
            if len(gene_name) > 0:
                gene_names.append(gene_name.values[0].split("_ENS")[0])
            else:
                gene_names.append(None)
        else:
            gene_names.append(None)
            top_sig.append(False)

    data_T["gene_name"] = gene_names
    data_T["top_sig"] = top_sig
    big_data_genes = data_T[data_T.gene_name.notnull()].copy()
    big_data_genes.set_index("gene_name", inplace=True)
    big_data_genes_filt = big_data_genes[big_data_genes.top_sig == True].copy()
    big_data_genes_filt.drop("top_sig", axis=1, inplace=True)
    return big_data_genes_filt


if __name__ == "__main__":
    deg_folder = sys.argv[1]
    top_n_sig = sys.argv[2].split(",")
    main(Path(deg_folder), top_n_sig)
