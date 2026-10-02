from datetime import datetime
import pandas as pd
from pathlib import Path
import json
import sys
from random import sample, shuffle


def main(prepare_folder):
    results_dir = Path(f"prompts/{str(datetime.now()).replace(' ','_')}")
    results_dir.mkdir(exist_ok=True)
    for file_name in prepare_folder.iterdir():
        if file_name.is_file():
            if file_name.suffix == ".json":
                print(file_name)
                with open(file_name, "r") as fp:
                    metadata = json.load(fp)

                train_out_fnames = []
                for train_annotated_filename in metadata["train_annotated_filenames"]:
                    model_data = pd.read_csv(train_annotated_filename)
                    out_fname = generate_prompts(
                        model_data, results_dir, train_annotated_filename
                    )
                    train_out_fnames.append(str(out_fname))

                holdout_out_fnames = []
                for holdout_annotated_filename in metadata[
                    "holdout_annotated_filenames"
                ]:
                    model_data = pd.read_csv(holdout_annotated_filename)
                    out_fname = generate_prompts(
                        model_data, results_dir, holdout_annotated_filename
                    )
                    holdout_out_fnames.append(str(out_fname))
                metadata["prompt_train_fnames"] = train_out_fnames
                metadata["prompt_holdout_fnames"] = holdout_out_fnames

                with open(results_dir / f"{file_name.stem}.json", "w") as fp:
                    json.dump(metadata, fp)


def generate_prompts(model_data, results_dir, train_annotated_filename):
    """
    early_mid_late_data = label_data[
        label_data.recurrence_time_sur.isin(["early", "mid", "late"])
    ]
    early_mid_late_panc_data = early_mid_late_data[early_mid_late_data.is_panc == True]
    print("early_mid_late_panc_data shape:", early_mid_late_panc_data.shape)
    """
    prompts = []
    acc_ids = []
    # df_sig_1k = model_data[model_data.top_1000_sig == True].copy()
    # df_sig_1k = df_sig_1k[df_sig_1k.gene_name.notnull()].copy()
    model_data.set_index("accession_id", inplace=True)
    label_cols = [
        "pfsmofromdx",
        "osmofromdx",
        "recurrence_time_sur",
    ]
    gene_data = model_data.drop(label_cols, axis=1)
    # big_data_genes = df_sig_1k.drop(columns=["top_1000_sig", "top_500_sig"])
    gene_count = 0
    for i in range(len(gene_data)):
        df_acc_id = gene_data.iloc[i].copy()
        df_acc_id.sort_values(ascending=False, inplace=True)
        gene_count = len(df_acc_id)

        gene_list = list(df_acc_id[:gene_count].index)
        cell_sentence = " ".join(sample(gene_list, gene_count))
        # cell_sentence = "MALAT1 TMSB4X B2M EEF1A1 H3F3B ACTB FTL RPL13 ..." # Truncated for example, use at least 200 genes for inference

        organism = "Homo sapiens"

        prompt = f"""The following is a list of {gene_count} gene names ordered by descending expression level in a {organism} cell. Your task is to give the cell type which this cell belongs to based on its gene expression.
        Cell sentence: {cell_sentence}.
        The cell type corresponding to these genes is:"""
        # print(f"prompt: {prompt}")
        prompts.append(prompt)
        acc_ids.append(df_acc_id.name)
    prompt_df = pd.DataFrame({"accession_id": acc_ids, "prompt": prompts})
    prompt_df_with_outcomes = prompt_df.merge(
        model_data.reset_index()[["accession_id"] + label_cols],
        on="accession_id",
        how="left",
    )
    in_fname = Path(train_annotated_filename)
    out_fname_stem = in_fname.name.split("_for_c2s_randomized")[0]
    out_fname = results_dir / f"prompts_{out_fname_stem}.csv"
    prompt_df_with_outcomes.to_csv(
        out_fname,
        index=False,
    )
    return out_fname


if __name__ == "__main__":
    prepare_folder = sys.argv[1]
    main(Path(prepare_folder))
