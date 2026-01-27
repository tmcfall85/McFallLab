from datetime import datetime
import pandas as pd
from random import sample

label_data = pd.read_csv("/mnt/c/Users/msochor/Downloads/big_data_with_labels.csv")
model_data = pd.read_csv("/mnt/c/Users/msochor/Downloads/big_data_for_c2s_modeling.csv")

early_mid_late_data = label_data[
    label_data.recurrence_time_sur.isin(["early", "mid", "late"])
]
early_mid_late_panc_data = early_mid_late_data[early_mid_late_data.is_panc == True]
print("early_mid_late_panc_data shape:", early_mid_late_panc_data.shape)
prompts = []
acc_ids = []
model_data = model_data[model_data.gene_name.notnull()].copy()
model_data.set_index("gene_name", inplace=True)
big_data_genes = model_data.drop(columns=["top_1000_sig", "top_500_sig"])
gene_count = 0
for c in big_data_genes.columns:
    if c in early_mid_late_panc_data.accession_id.values:
        # print(f"Acc id: {c}")
        df_acc_id = big_data_genes[c].copy()
        df_acc_id.sort_values(ascending=False, inplace=True)

        gene_count = 500

        gene_list = list(df_acc_id[:gene_count].index)
        cell_sentence = " ".join(sample(gene_list, gene_count))
        # cell_sentence = "MALAT1 TMSB4X B2M EEF1A1 H3F3B ACTB FTL RPL13 ..." # Truncated for example, use at least 200 genes for inference

        organism = "Homo sapiens"

        prompt = f"""The following is a list of {gene_count} gene names ordered by descending expression level in a {organism} cell. Your task is to give the cell type which this cell belongs to based on its gene expression.
        Cell sentence: {cell_sentence}.
        The cell type corresponding to these genes is:"""
        # print(f"prompt: {prompt}")
        prompts.append(prompt)
        acc_ids.append(c)
prompt_df = pd.DataFrame({"accession_id": acc_ids, "prompt": prompts})
out_fname = f"prompts_{gene_count}_highest_tpm_randomized_{str(datetime.now()).replace(' ','_')}"
prompt_df.to_csv(
    f"{out_fname}.csv",
    index=False,
)
print(len(prompt_df))
print(f"filename: {out_fname}")
