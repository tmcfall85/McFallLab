import shutil
from pathlib import Path, PosixPath

import pandas as pd


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
            bam_file = acc_rna_folder / f"{accession_number}_alig_csort.bam"
            if bam_file.is_file():
                bam_file_exists = True
                dest_file = scratch_dir / f"{accession_number}_alig_csort.bam"
                shutil.copy(str(bam_file), str(dest_file))
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
            tempus_dirs.append(tempus_dir)

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
    move_rcc_files_to_scratch(
        ["TL-19-695788", "TL-19-4A2D28"],
        Path("/group/dseo"),
        ["seo_tempus_first_seq_data", "seo_tempus_second_seq_data"],
        Path("/scratch/g/dseo/rna_seq_temp/"),
    )
