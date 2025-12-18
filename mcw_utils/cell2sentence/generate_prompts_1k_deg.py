from datetime import datetime
import pandas as pd

label_data = pd.read_csv("/mnt/c/Users/msochor/Downloads/big_data_with_labels.csv")
model_data = pd.read_csv("/mnt/c/Users/msochor/Downloads/big_data_for_c2s_modeling.csv")

early_mid_late_data = label_data[
    label_data.recurrence_time_sur.isin(["early", "mid", "late"])
]
early_mid_late_panc_data = early_mid_late_data[early_mid_late_data.is_panc == True]
print("early_mid_late_panc_data shape:", early_mid_late_panc_data.shape)
prompts = []
acc_ids = []
df_sig_1k = model_data[model_data.top_1000_sig == True].copy()
df_sig_1k = df_sig_1k[df_sig_1k.gene_name.notnull()].copy()
df_sig_1k.set_index("gene_name", inplace=True)
big_data_genes = df_sig_1k.drop(columns=["top_1000_sig"])
gene_count = 0
for c in big_data_genes.columns:
    if c in early_mid_late_panc_data.accession_id.values:
        # print(f"Acc id: {c}")
        df_acc_id = big_data_genes[c].copy()
        df_acc_id.sort_values(ascending=False, inplace=True)

        gene_count = len(big_data_genes)

        cell_sentence = " ".join(df_acc_id[:gene_count].index)
        # cell_sentence = "MALAT1 TMSB4X B2M EEF1A1 H3F3B ACTB FTL RPL13 ..." # Truncated for example, use at least 200 genes for inference

        organism = "Homo sapiens"

        prompt = f"""The following is a list of {gene_count} gene names ordered by descending expression level in a {organism} cell. Your task is to give the cell type which this cell belongs to based on its gene expression.
        Cell sentence: {cell_sentence}.
        The cell type corresponding to these genes is:"""
        # print(f"prompt: {prompt}")
        prompts.append(prompt)
        acc_ids.append(c)
prompt_df = pd.DataFrame({"accession_id": acc_ids, "prompt": prompts})
out_fname = f"prompts_{gene_count}_deg_{str(datetime.now()).replace(' ','_')}"
prompt_df.to_csv(
    f"{out_fname}.csv",
    index=False,
)
print(len(prompt_df))
print(f"filename: {out_fname}")
