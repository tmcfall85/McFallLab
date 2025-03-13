import shutil
from pathlib import Path, PosixPath

import pandas as pd
import sys


def move_rcc_file_to_scratch(
    accession_number: str, group_dir: PosixPath, scratch_dir: PosixPath
):

    acc_folder = group_dir / accession_number
    print(acc_folder)
    print(acc_folder.is_dir())
    if acc_folder.is_dir():
        acc_number_exists = True
        acc_rna_folder = acc_folder / "RNA"
        if acc_rna_folder.is_dir():
            rna_folder_exists = True
            for item in acc_rna_folder.iterdir():
                if item.suffix == '.bam' and item.is_file():
                    bam_file_exists = True
                    dest_file = scratch_dir / item.name
                    shutil.copy(str(item), str(dest_file))
            else:
                bam_file_exists = False
        else:
            rna_folder_exists = False
            bam_file_exists = False
    else:
        acc_number_exists = False
        rna_folder_exists = False
        bam_file_exists = False

    return acc_number_exists, rna_folder_exists, bam_file_exists


def move_rcc_files_to_scratch(
    accession_numbers: list[str],
    group_dir: PosixPath,
    group_subdirs: list[str],
    scratch_dir: PosixPath,
):
    acc_numbers_exists = []
    rna_folders_exists = []
    bam_files_exists = []
    tempus_dirs = []
    for accession_number in accession_numbers:
        for group_subdir in group_subdirs:
            tempus_dir = group_dir / group_subdir
            acc_number_exists, rna_folder_exists, bam_file_exists = (
                move_rcc_file_to_scratch(accession_number, tempus_dir, scratch_dir)
            )
            acc_numbers_exists.append(acc_number_exists)
            rna_folders_exists.append(rna_folder_exists)
            bam_files_exists.append(bam_file_exists)
            tempus_dirs.append(tempus_dir / accession_number)

    move_report = pd.DataFrame(
        {
            "tempus_dir": tempus_dirs,
            "acc_number_exists": acc_numbers_exists,
            "rna_folder_exists": rna_folders_exists,
            "bam_file_exists": bam_files_exists,
        }
    )
    move_report.to_csv(scratch_dir / "move_report.csv", index=False)


if __name__ == "__main__":
    user = sys.argv[1]
    acc_file = pd.read_csv(sys.argv[2])
    out_dir = sys.argv[3]
    sub_dirs = sys.argv[4:]

    move_rcc_files_to_scratch(
        acc_file.accession_id.values,
        Path(f"/group/{user}"),
        sub_dirs,
        Path(f"/scratch/g/{user}/rna_seq_temp/"),
    )
