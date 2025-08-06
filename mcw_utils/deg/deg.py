import pandas as pd
import sys


def main(deg_file, dta_file, merge_file):
    deg_df = pd.read_csv(deg_file)
    dta_df = pd.read_stata(dta_file)
    merge_df = pd.read_csv(merge_file)

    # Perform analysis or processing
    merge_has_mrn = merge_df[merge_df.mrn.notnull()]
    unique_acc_mrns = (
        merge_has_mrn.groupby(["accession_id", "mrn", "emrn"])
        .emr_id.count()
        .reset_index()
    )

    deg_with_mrn = deg_df.merge(unique_acc_mrns, on="accession_id", how="left")
    deg_with_mrn = deg_with_mrn[deg_with_mrn.mrn.notnull()]
    deg_with_pfs = deg_with_mrn.merge(
        dta_df[["mrn", "recurrence_time_sur"]], on="mrn", how="left"
    )
    print("Distribution of recurrence_time_sur:")
    print(deg_with_pfs.groupby("recurrence_time_sur").accession_id.count())
    return


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python deg.py <deg_file> <dta_file> <merge_file>")
        sys.exit(1)
    deg_file = sys.argv[1]
    dta_file = sys.argv[2]
    merge_file = sys.argv[3]
    obj = main(deg_file, dta_file, merge_file)
