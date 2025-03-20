import sys

import numpy as np
import pandas as pd


def any_rna_raw(df):
    return np.any(df["Has RNA Folder"] == "Yes")


def any_dna_raw(df):
    return np.any(df["Has DNA Folder"] == "Yes")


def any_rna_json(df):
    return np.any(df.report_type == "RNA")


def any_dna_json(df):
    return np.any(df.report_type == "DNA")


def label_by_bmi(pdac_sotb, patients):
    inner_patients = patients.copy()
    inner_patients["mrn"] = None
    for i in range(len(inner_patients)):
        filt = pdac_sotb[pdac_sotb.bmi == inner_patients.iloc[i]["BMI at diagnosis"]]
        filt = filt[filt.bmi_restaging == inner_patients.iloc[i]["BMI at pre-op"]]
        if len(filt) == 1:
            inner_patients.loc[i, "mrn"] = filt.iloc[0].mrn
        else:
            if np.isnan(inner_patients.iloc[i]["Tissue Bank Number"]) == False:
                print(f"No unique BMI label found for {inner_patients.iloc[i]}")
    return inner_patients


def label(merged_tempus_fname: str, patients_fname: str, pdac_sotb_fname: str):

    merged_tempus = pd.read_csv(merged_tempus_fname)
    patients = pd.read_excel(patients_fname)
    pdac_sotb = pd.read_stata(pdac_sotb_fname)

    patients_labeled = label_by_bmi(pdac_sotb, patients)

    control = patients_labeled.loc[:19, :]
    test = patients_labeled.loc[22:, :]

    control_merged = control.merge(merged_tempus, on="mrn", how="left")
    test_merged = test.merge(merged_tempus, on="mrn", how="left")

    test_dna_raw = (
        test_merged.groupby("mrn")
        .apply(any_dna_raw)
        .reset_index()
        .rename({0: "DNA_data_on_rcc"}, axis=1)
    )
    test_rna_raw = (
        test_merged.groupby("mrn")
        .apply(any_rna_raw)
        .reset_index()
        .rename({0: "RNA_data_on_rcc"}, axis=1)
    )
    control_dna_raw = (
        control_merged.groupby("mrn")
        .apply(any_dna_raw)
        .reset_index()
        .rename({0: "DNA_data_on_rcc"}, axis=1)
    )
    control_rna_raw = (
        control_merged.groupby("mrn")
        .apply(any_rna_raw)
        .reset_index()
        .rename({0: "RNA_data_on_rcc"}, axis=1)
    )

    test_dna_json = (
        test_merged.groupby("mrn")
        .apply(any_dna_json)
        .reset_index()
        .rename({0: "DNA_json"}, axis=1)
    )
    test_rna_json = (
        test_merged.groupby("mrn")
        .apply(any_rna_json)
        .reset_index()
        .rename({0: "RNA_json"}, axis=1)
    )
    control_dna_json = (
        control_merged.groupby("mrn")
        .apply(any_dna_json)
        .reset_index()
        .rename({0: "DNA_json"}, axis=1)
    )
    control_rna_json = (
        control_merged.groupby("mrn")
        .apply(any_rna_json)
        .reset_index()
        .rename({0: "RNA_json"}, axis=1)
    )

    test_readout = (
        test.merge(test_dna_raw, on="mrn")
        .merge(test_rna_raw, on="mrn")
        .merge(test_dna_json, on="mrn")
        .merge(test_rna_json, on="mrn")
    )
    control_readout = (
        control.merge(control_dna_raw, on="mrn")
        .merge(control_rna_raw, on="mrn")
        .merge(control_dna_json, on="mrn")
        .merge(control_rna_json, on="mrn")
    )

    test_readout.to_csv("test_data_rna_dna_missingness.csv", index=False)
    control_readout.to_csv("control_data_rna_dna_missingness.csv", index=False)
    control_merged.to_csv("control_merged_to_tempus_data.csv", index=False)
    test_merged.to_csv("test_merged_to_tempus_data.csv", index=False)


if __name__ == "__main__":
    label(sys.argv[1], sys.argv[2], sys.argv[3])
