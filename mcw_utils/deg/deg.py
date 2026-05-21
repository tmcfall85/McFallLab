import pandas as pd
import sys
import requests
import time
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from datetime import date, datetime
from tqdm import tqdm
import numpy as np
import pickle
from pathlib import Path
import json

# Create a single directory

from sklearn.model_selection import StratifiedShuffleSplit, train_test_split
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
    tempus_merge = pd.read_csv(merge_file)
    print("Input files read successfully.")

    merge_columns = ['patient_mrn', 'accession_id', 'sample_site']
    print("Merging dataframes...")
    tempus_merge['sample_site'] = tempus_merge.sample_site.fillna('')
    print(f'Merge file length: {len(tempus_merge)}')
    tempus_merge_unique = tempus_merge.groupby(merge_columns).count().reset_index()
    
    print(f'Length after removing duplicates: {len(tempus_merge_unique)}')
    tempus_deg_df = deg_df.merge(tempus_merge_unique[merge_columns], on='accession_id', how='inner')
    print(f"Merged deg_df (len: {len(deg_df)}) to tempus_merge_unique (len: {len(tempus_merge_unique)}): result len {len(tempus_deg_df)}")
    pdac_tempus_deg_df = tempus_deg_df[tempus_deg_df.sample_site.str.startswith('Pancreas')].copy()
    print(f'Length after PDAC filter: {len(pdac_tempus_deg_df)}')
    pdac_tempus_deg_df.drop('sample_site', axis=1, inplace=True)
    #merge_has_mrn = merge_df[merge_df.mrn.notnull()]
    #merge_has_mrn = merge_has_mrn[merge_has_mrn.report_type == "RNA"]
    #unique_acc_mrns = (
    #    merge_has_mrn.groupby(["accession_id", "mrn", "emrn", "specimen_sample_site"])
    #    .emr_id.count()
    #    .reset_index()
    #)
    #unique_acc_mrns.drop("emr_id", axis=1, inplace=True)
    dta_merge_columns = ['mrn', 'pfsmofromdx', 'osmofromdx', 'recurrence_time_sur']
    #deg_with_mrn = deg_df.merge(unique_acc_mrns, on="accession_id", how="right")
    # deg_with_mrn = deg_with_mrn[deg_with_mrn.mrn.notnull()]
    deg_with_pfs = pdac_tempus_deg_df.merge(
        dta_df[dta_merge_columns], left_on="patient_mrn", right_on="mrn", how="inner"
    )
    print(f"Merged pdac_tempus_deg_df (len: {len(pdac_tempus_deg_df)}) to dta_df (len: {len(dta_df)}): result len {len(deg_with_pfs)}")
    
    deg_with_pfs.drop(columns=["patient_mrn", "mrn"], inplace=True)
    deg_with_pfs_eml = deg_with_pfs[
        deg_with_pfs.recurrence_time_sur.isin(["early", "mid", "late"])
    ].copy()
    print(f'Filtered to just early/mid/late. Len before {len(deg_with_pfs)}, len after {len(deg_with_pfs_eml)}')
   
    print("Distribution of recurrence_time_sur after merges:")
    print(deg_with_pfs_eml.groupby(["recurrence_time_sur"]).accession_id.count())
    results_dir = Path(f"results/{str(datetime.now()).replace(' ','_')}")
    results_dir.mkdir(exist_ok=True)
    deg_with_pfs_eml.to_csv(results_dir / "deg_with_pfs.csv", index=False)
    metadata = {
        'deg_file': deg_file, 
        'dta_file': dta_file, 
        'merge_file': merge_file
    }
    with open(results_dir / 'metadata.json', 'w') as fp: 
        json.dump(metadata, fp)
    return deg_with_pfs_eml, results_dir



def kfold_run_deseq2(
    deg_with_pfs, results_dir, holdout_fraction, a_label="early", b_label="late", n_splits=5
):

    
    # train_ds, val_ds = random_split(dataset, [n_train, n_val])

    # Assuming you have your features X and target y as numpy arrays or pandas DataFrames/Series
    # X and y must have the same number of samples (e.g., n_samples = 200)
    # Replace n_samples with the actual number of samples in your dataset

    splits = {}
    if holdout_fraction > 0:
        X_train, X_holdout, y_train, y_holdout = train_test_split(deg_with_pfs, deg_with_pfs.recurrence_time_sur,
                                                        stratify=deg_with_pfs.recurrence_time_sur, 
                                                        test_size=holdout_fraction)
        X_holdout.to_csv(results_dir / 'holdout.csv', index=False)
        X_train.to_csv(results_dir / 'train.csv', index=False)
        print(f'Holdout saved with len {len(X_holdout)}')
        print(f'Train len {len(X_train)}')
    else:
        X_train = deg_with_pfs
        y_train = deg_with_pfs.recurrence_time_sur

    n_val_max = np.ceil(len(X_train) * (1 / n_splits))
    n_val_min = np.floor(len(X_train) * (1 / n_splits))
    print(n_val_min, n_val_max)
    if (len(X_train) - n_val_max) % 2 == 0:
        n_val = int(n_val_max)
    else:
        n_val = int(n_val_min)
    print(len(X_train) - n_val, n_val)

    n_train = len(X_train) - n_val
    
    ss = StratifiedShuffleSplit(n_splits=n_splits, test_size=n_val)

    for i, (train_index, val_index) in enumerate(
        ss.split(X_train, y_train)
    ):
        print(f"Fold {i+1}:")
        print(f"  Train set size: {len(train_index)}")
        print(f"  Validation set size: {len(val_index)}")
        splits[i] = (train_index, val_index)

        deg_with_pfs_train = deg_with_pfs.iloc[train_index]
        run_deseq2(
            deg_with_pfs_train,
            results_dir,
            f"deg_results_kfold_{i}",
            a_label=a_label,
            b_label=b_label,
        )

    with open(results_dir / f"fold_indices.pkl", "wb") as fp:
        pickle.dump(splits, fp)


def run_deseq2(deg_with_pfs, results_dir, outfile_name, a_label="early", b_label="late"):
    print(f"Running DESeq2 analysis on conditions: {a_label} vs {b_label}...")
    deg_with_pfs_ab = deg_with_pfs[
        deg_with_pfs.recurrence_time_sur.isin([a_label, b_label])
    ].copy()

    deg_with_pfs_ab.set_index("accession_id", inplace=True)
    metadata = deg_with_pfs_ab[["recurrence_time_sur", "pfsmofromdx"]].copy()
    
    a_condition = lambda row: condition(row, a_label)
    metadata["condition"] = metadata.apply(a_condition, axis=1)
    metadata.drop(columns=["pfsmofromdx"], inplace=True)


    deg_with_pfs_ab.drop(columns=["recurrence_time_sur", "pfsmofromdx", "osmofromdx"], inplace=True)
    deg_with_pfs_ab_int = deg_with_pfs_ab.copy()
    print("Converting TPMs to integers for deseq...")
    for col in tqdm(deg_with_pfs_ab_int.columns):
        deg_with_pfs_ab_int[col] = integerize(deg_with_pfs_ab_int[col])

    print("Filtering genes with less than 1000 counts across samples...")
    genes_to_keep = deg_with_pfs_ab_int.columns[deg_with_pfs_ab_int.sum(axis=0) >= 1000]
    print("Genes before filtering:", deg_with_pfs_ab_int.shape[1])
    deg_with_pfs_ab_int = deg_with_pfs_ab_int[genes_to_keep]
    print("Genes after filtering:", deg_with_pfs_ab_int.shape[1])
    print(f"count shape: {deg_with_pfs_ab_int.shape}")
    print(f'metadata shape: {metadata.shape}')
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
    outfile = results_dir / f"{outfile_name}.csv"
    results_with_metadata.to_csv(outfile)
    print(f"Results saved to {outfile}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            "Usage: python deg.py <deg_file> <dta_file> <merge_file> <k_fold> <holdout_fraction>"
        )
        sys.exit(1)
    deg_file = sys.argv[1]
    dta_file = sys.argv[2]
    merge_file = sys.argv[3]
    use_kfold = sys.argv[4] if len(sys.argv) > 4 else "false"
    holdout_fraction = float(sys.argv[5]) if len(sys.argv) > 5 else 0

    merged, results_dir = read_and_merge(deg_file, dta_file, merge_file)
    if use_kfold.lower() == "true":
        kfold_run_deseq2(merged, results_dir, holdout_fraction)
    else:
        run_deseq2(merged, results_dir, "deg_results")
