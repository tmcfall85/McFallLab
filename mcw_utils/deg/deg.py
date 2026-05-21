import pandas as pd
import sys
import requests
import time
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from datetime import date
from tqdm import tqdm
import numpy as np
import pickle

from sklearn.model_selection import StratifiedShuffleSplit

# Reproducibility
RNG_SEED = 42
np.random.seed(RNG_SEED)


def integerize(col):
    x = 100 * col
    return x.astype(int)


def condition(row, a_label):
    if row.recurrence_time_sur == a_label:
        return "A"
    else:
        return "B"


def get_gene_metadata(gene_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=%28{gene_id}%29&fields=protein_name%2C%20gene_names%2C%20cc_function&size=5"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()
    # Extract the protein name from the first result, if available
    try:
        return (
            data["results"][0]["proteinDescription"]["recommendedName"]["fullName"][
                "value"
            ],
            data["results"][0]["genes"][0]["geneName"]["value"],
            data["results"][0]["comments"][0]["texts"][0]["value"],
        )
    except (KeyError, IndexError):
        return None, None, None


def read_and_merge(deg_file, dta_file, merge_file):
    # Read input files
    print("Reading input files...")
    deg_df = pd.read_csv(deg_file)
    dta_df = pd.read_stata(dta_file)
    merge_df = pd.read_csv(merge_file)
    print("Input files read successfully.")

    print("Merging dataframes...")
    merge_has_mrn = merge_df[merge_df.mrn.notnull()]
    merge_has_mrn = merge_has_mrn[merge_has_mrn.report_type == "RNA"]
    unique_acc_mrns = (
        merge_has_mrn.groupby(["accession_id", "mrn", "emrn", "specimen_sample_site"])
        .emr_id.count()
        .reset_index()
    )
    unique_acc_mrns.drop("emr_id", axis=1, inplace=True)

    deg_with_mrn = deg_df.merge(unique_acc_mrns, on="accession_id", how="right")
    # deg_with_mrn = deg_with_mrn[deg_with_mrn.mrn.notnull()]
    deg_with_pfs = deg_with_mrn.merge(
        dta_df[["mrn", "recurrence_time_sur"]], on="mrn", how="left"
    )
    deg_with_pfs.drop(columns=["mrn", "emrn"], inplace=True)
    deg_with_pfs.specimen_sample_site.fillna("unknown", inplace=True)
    deg_with_pfs["is_panc"] = deg_with_pfs.apply(is_panc, axis=1)
    deg_with_pfs_panc = deg_with_pfs[deg_with_pfs.is_panc == True].copy()
    deg_with_pfs_panc_eml = deg_with_pfs_panc[
        deg_with_pfs_panc.recurrence_time_sur.isin(["early", "mid", "late"])
    ].copy()
    deg_with_pfs_panc_eml.drop(
        columns=["is_panc", "specimen_sample_site"], inplace=True
    )
    print("Distribution of recurrence_time_sur after merges:")
    print(deg_with_pfs_panc_eml.groupby(["recurrence_time_sur"]).accession_id.count())
    deg_with_pfs_panc_eml.to_csv("deg_with_pfs.csv")
    return deg_with_pfs_panc_eml


def is_panc(row):
    return row.specimen_sample_site.startswith("Panc")


def kfold_run_deseq2(
    deg_with_pfs, out_file, a_label="early", b_label="late", n_splits=5
):

    n_val_max = np.ceil(len(deg_with_pfs) * (1 / n_splits))
    n_val_min = np.floor(len(deg_with_pfs) * (1 / n_splits))
    print(n_val_min, n_val_max)
    if (len(deg_with_pfs) - n_val_max) % 2 == 0:
        n_val = int(n_val_max)
    else:
        n_val = int(n_val_min)
    print(len(deg_with_pfs) - n_val, n_val)

    n_train = len(deg_with_pfs) - n_val
    # train_ds, val_ds = random_split(dataset, [n_train, n_val])

    # Assuming you have your features X and target y as numpy arrays or pandas DataFrames/Series
    # X and y must have the same number of samples (e.g., n_samples = 200)
    # Replace n_samples with the actual number of samples in your dataset

    ss = StratifiedShuffleSplit(n_splits=n_splits, test_size=n_val)
    splits = {}
    for i, (train_index, val_index) in enumerate(
        ss.split(deg_with_pfs, deg_with_pfs.recurrence_time_sur)
    ):
        print(f"Fold {i+1}:")
        print(f"  Train set size: {len(train_index)}")
        print(f"  Validation set size: {len(val_index)}")
        splits[i] = (train_index, val_index)

        deg_with_pfs_train = deg_with_pfs.iloc[train_index]
        run_deseq2(
            deg_with_pfs_train,
            f"{out_file}_kfold_{i}",
            a_label=a_label,
            b_label=b_label,
        )

    with open(f"{out_file}_fold_indices.pkl", "wb") as fp:
        pickle.dump(splits, fp)


def run_deseq2(deg_with_pfs, out_file, a_label="early", b_label="late"):
    print(f"Running DESeq2 analysis on conditions: {a_label} vs {b_label}...")
    deg_with_pfs_ab = deg_with_pfs[
        deg_with_pfs.recurrence_time_sur.isin([a_label, b_label])
    ].copy()

    deg_with_pfs_ab.set_index("accession_id", inplace=True)

    metadata = deg_with_pfs_ab[["recurrence_time_sur", "emr_id"]].copy()
    a_condition = lambda row: condition(row, a_label)
    metadata["condition"] = metadata.apply(a_condition, axis=1)
    metadata.drop(columns=["emr_id"], inplace=True)

    deg_with_pfs_ab.drop(columns=["recurrence_time_sur", "emr_id"], inplace=True)
    deg_with_pfs_ab_int = deg_with_pfs_ab.copy()
    print("Converting TPMs to integers for deseq...")
    for col in tqdm(deg_with_pfs_ab_int.columns):
        deg_with_pfs_ab_int[col] = integerize(deg_with_pfs_ab_int[col])

    print("Filtering genes with less than 1000 counts across samples...")
    genes_to_keep = deg_with_pfs_ab_int.columns[deg_with_pfs_ab_int.sum(axis=0) >= 1000]
    print("Genes before filtering:", deg_with_pfs_ab_int.shape[1])
    deg_with_pfs_ab_int = deg_with_pfs_ab_int[genes_to_keep]
    print("Genes after filtering:", deg_with_pfs_ab_int.shape[1])
    inference = DefaultInference(n_cpus=8)
    dds = DeseqDataSet(
        counts=deg_with_pfs_ab_int,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
    )
    print("Running DESeq2...")
    dds.deseq2()

    ds = DeseqStats(dds, contrast=["condition", "B", "A"], inference=inference)

    ds.summary()
    results_df = ds.results_df.copy()
    results_df["baseMean"] = results_df.baseMean / 100
    print("Labeling top results...")
    results_df["logfold_sort"] = 1 / results_df.log2FoldChange.abs()
    full_names = []
    gene_names = []
    comments = []
    indices = []
    for index in tqdm(results_df[results_df.padj < 0.05].index):
        time.sleep(1)
        full_name, gene_name, comment = get_gene_metadata(index.split(".")[0])
        full_names.append(full_name)
        gene_names.append(gene_name)
        comments.append(comment)
        indices.append(index)
    gene_metadata = pd.DataFrame(
        {"full_name": full_names, "gene_name": gene_names, "comment": comments},
        index=indices,
    )
    results_with_metadata = pd.merge(
        results_df, gene_metadata, left_index=True, right_index=True, how="left"
    )
    results_with_metadata.sort_values(
        by=["padj", "logfold_sort"], ascending=True, inplace=True
    )
    results_with_metadata.drop(columns=["logfold_sort"], inplace=True)
    results_with_metadata.fillna("", inplace=True)
    out_file_today = f"{out_file}_{date.today()}.csv"
    results_with_metadata.to_csv(out_file_today)
    print(f"Results saved to {out_file_today}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            "Usage: python deg.py <deg_file> <dta_file> <merge_file> <k_fold> <out_file>"
        )
        sys.exit(1)
    deg_file = sys.argv[1]
    dta_file = sys.argv[2]
    merge_file = sys.argv[3]
    use_kfold = sys.argv[4] if len(sys.argv) > 4 else "false"
    out_file = sys.argv[5] if len(sys.argv) > 5 else "deg_results_with_metadata"

    merged = read_and_merge(deg_file, dta_file, merge_file)
    if use_kfold.lower() == "true":
        kfold_run_deseq2(merged, out_file)
    else:
        run_deseq2(merged, out_file)
