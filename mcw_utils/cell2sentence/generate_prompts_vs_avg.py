from datetime import datetime
import pandas as pd

label_data = pd.read_csv("/mnt/c/Users/msochor/Downloads/big_data_with_labels.csv")
model_data = pd.read_csv("/mnt/c/Users/msochor/Downloads/big_data_for_c2s_modeling.csv")
avg_label = "late"

avg_mid_late_data = label_data[
    label_data.recurrence_time_sur.isin(["early", "mid", "late"])
]
avg_mid_late_panc_data = avg_mid_late_data[avg_mid_late_data.is_panc == True]
print(f"{avg_label}_mid_late_panc_data shape:", avg_mid_late_panc_data.shape)

avg_data = label_data[label_data.recurrence_time_sur.isin([avg_label])]
avg_panc_data = avg_data[avg_data.is_panc == True]
print(f"{avg_label}_mid_late_panc_data shape:", avg_panc_data.shape)
prompts = []
acc_ids = []

df_sig_500 = model_data[model_data.top_500_sig == True].copy()
df_sig_500 = df_sig_500[df_sig_500.gene_name.notnull()].copy()
df_sig_500.set_index("gene_name", inplace=True)
big_data_genes = df_sig_500.drop(columns=["top_500_sig"])
gene_count = 0
drop_rows = []
for c in big_data_genes.columns:
    if c not in avg_panc_data.accession_id.values:
        drop_rows.append(c)
avg_big_data_genes = big_data_genes.drop(columns=drop_rows)
sorted_avg_genes = avg_big_data_genes.mean(axis=1).sort_values(ascending=False)
avg_cell_sentence = " ".join(sorted_avg_genes.index)
gene_count = 0
for c in big_data_genes.columns:
    if c in avg_mid_late_panc_data.accession_id.values:
        # print(f"Acc id: {c}")
        df_acc_id = big_data_genes[c].copy()
        df_acc_id.sort_values(ascending=False, inplace=True)

        gene_count = len(big_data_genes)

        cell_sentence = " ".join(df_acc_id[:gene_count].index)
        # cell_sentence = "MALAT1 TMSB4X B2M EEF1A1 H3F3B ACTB FTL RPL13 ..." # Truncated for example, use at least 200 genes for inference

        organism = "Homo sapiens"

        prompt = f"""Given data derived from pancreatic cells from Homo Sapiens and the following cell sentence of {gene_count} expressed genes representing a cell's basal state and the cell sentence of {gene_count} expressed genes representing a cell's perturbed state, predict the perturbation.
        Control cell sentence: {avg_cell_sentence}.
        
        Perturbed cell sentence: {cell_sentence}.
        
        The pertubation is: """
        # print(f"prompt: {prompt}")
        prompts.append(prompt)
        acc_ids.append(c)
prompt_df = pd.DataFrame({"accession_id": acc_ids, "prompt": prompts})
out_fname = f"prompts_vs_avg_{gene_count}_{str(datetime.now()).replace(' ','_')}"
prompt_df.to_csv(
    f"{out_fname}.csv",
    index=False,
)
print(len(prompt_df))
print(f"filename: {out_fname}")
