import shutil
from pathlib import Path, PosixPath

import pandas as pd
import sys


def move_rcc_file_to_scratch(
    accession_number: str, group_dir: PosixPath, scratch_dir: PosixPath
):

    acc_folder = group_dir
    if acc_folder.is_dir():
        acc_number_exists = True
        acc_star_folder = acc_folder / "star_output"
        if acc_star_folder.is_dir():
            rna_folder_exists = True
            for item in acc_star_folder.iterdir():
                if (
                    item.suffix == ".bam"
                    and item.is_file()
                    and item.name.find("sortedByCoord") > -1
                ):
                    bam_file_exists = True
                    dest_name = f"{accession_number}_alig.bam"
                    dest_file = scratch_dir / dest_name
                    print(f"copying {item} to {dest_file}")
                    shutil.copy(str(item), str(dest_file))
            else:
                bam_file_exists = False
        else:
            rna_folder_exists = False
            bam_file_exists = False
        acc_rsem_folder = acc_folder / "rsem_output"
        if acc_rsem_folder.is_dir():
            for item in acc_rsem_folder.iterdir():
                if (
                    item.suffix == ".results"
                    and item.is_file()
                    and item.name.find("genes") > -1
                ):
                    rsem_file_exists = True
                    dest_name = f"{accession_number}_rsem.results"
                    dest_file = scratch_dir / dest_name
                    print(f"copying {item} to {dest_file}")
                    shutil.copy(str(item), str(dest_file))
            else:
                rsem_file_exists = False
        else:
            rsem_file_exists = False
    else:
        acc_number_exists = False
        rna_folder_exists = False
        bam_file_exists = False
        rsem_file_exists = False

    return acc_number_exists, rna_folder_exists, bam_file_exists, rsem_file_exists


def move_rcc_files_to_scratch(
    group_dir: PosixPath,
    scratch_dir: PosixPath,
):
    if scratch_dir.is_dir() == False:
        scratch_dir.mkdir()

    acc_numbers_exists = []
    rna_folders_exists = []
    bam_files_exists = []
    rsem_files_exisxts = []
    tempus_dirs = []
    for tempus_dir in group_dir.iterdir():
        accession_number = tempus_dir.stem
        if tempus_dir.is_dir():
            (
                acc_number_exists,
                rna_folder_exists,
                bam_file_exists,
                rsem_file_exists,
            ) = move_rcc_file_to_scratch(accession_number, tempus_dir, scratch_dir)
            acc_numbers_exists.append(acc_number_exists)
            rna_folders_exists.append(rna_folder_exists)
            bam_files_exists.append(bam_file_exists)
            rsem_files_exisxts.append(rsem_file_exists)
            tempus_dirs.append(tempus_dir)

    move_report = pd.DataFrame(
        {
            "tempus_dir": tempus_dirs,
            "acc_number_exists": acc_numbers_exists,
            "rna_folder_exists": rna_folders_exists,
            "bam_file_exists": bam_files_exists,
            "rsem_file_exists": rsem_files_exisxts,
        }
    )
    move_report.to_csv(scratch_dir / "move_report.csv", index=False)


if __name__ == "__main__":
    user = sys.argv[1]
    in_dir = sys.argv[2]
    out_dir = sys.argv[3]

    move_rcc_files_to_scratch(
        Path(f"/group/{user}/{in_dir}"),
        Path(f"/scratch/g/{user}/{out_dir}/"),
    )
