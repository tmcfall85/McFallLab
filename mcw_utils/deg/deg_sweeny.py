import pandas as pd
import sys
import requests
import time
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from datetime import date
from tqdm import tqdm


def integerize(col):
    x = 100 * col
    return x.astype(int)


def condition(row, a_label):
    if row.group == a_label:
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


def read_and_merge(deg_file, experiment_file):
    deg_df = pd.read_csv(deg_file)
    deg_df.rename(columns={"Unnamed: 0": "gene_id"}, inplace=True)
    deg_df.set_index("gene_id", inplace=True)
    experiment_df = pd.read_csv(experiment_file)

    deg_df = deg_df.T.reset_index()
    deg_df.rename(columns={"index": "experiment"}, inplace=True)
    deg_df["experiment"] = pd.to_numeric(deg_df["experiment"], errors="coerce").astype(
        "Int64"
    )
    deg_with_experiment = deg_df.merge(experiment_df, on="experiment", how="left")
    print(deg_with_experiment.head())
    return deg_with_experiment


def run_deseq2(deg_with_pfs, out_file, a_label="GFP", b_label="NOX5"):
    print(f"Running DESeq2 analysis on conditions: {a_label} vs {b_label}...")
    deg_with_pfs_ab = deg_with_pfs[deg_with_pfs.group.isin([a_label, b_label])].copy()
    metadata = deg_with_pfs_ab[["group", "experiment"]].copy()
    deg_with_pfs_ab.set_index("experiment", inplace=True)

    metadata = deg_with_pfs_ab[["group"]].copy()
    a_condition = lambda row: condition(row, a_label)
    metadata["condition"] = metadata.apply(a_condition, axis=1)
    print(metadata.head())
    metadata.drop(columns=["group"], inplace=True)

    deg_with_pfs_ab.drop(columns=["group"], inplace=True)
    deg_with_pfs_ab_int = deg_with_pfs_ab.copy()
    print("Converting TPMs to integers for deseq...")
    deg_with_pfs_ab_int.fillna(0, inplace=True)
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
    deg_to_merge = deg_with_pfs_ab_int.T
    results_with_metadata = results_with_metadata.merge(
        deg_to_merge, left_index=True, right_index=True, how="left"
    )
    out_file_today = f"{out_file}_{date.today()}.csv"
    results_with_metadata.to_csv(out_file_today)
    print(f"Results saved to {out_file_today}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(
            "Usage: python deg.py <deg_file> <experiment_file> <a_label> <b_label> <out_file>"
        )
        sys.exit(1)
    deg_file = sys.argv[1]
    experiment_file = sys.argv[2]
    a_label = sys.argv[3]
    b_label = sys.argv[4]
    out_file = sys.argv[5] if len(sys.argv) > 5 else "deg_results_with_metadata"
    merged = read_and_merge(deg_file, experiment_file)
    run_deseq2(merged, out_file, a_label=a_label, b_label=b_label)
