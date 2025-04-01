import shutil
from pathlib import Path, PosixPath

import pandas as pd
import sys
import tarfile


def extract_tar_gz(file_path, extract_path):
    try:
        with tarfile.open(file_path, "r:gz") as tar:
            tar.extractall(extract_path)
        print(f"Successfully extracted '{file_path}' to '{extract_path}'")
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except tarfile.ReadError:
        print(f"Error: '{file_path}' is not a valid tar.gz file.")
    except Exception as e:
        print(f"An error occurred: {e}")


def move_rcc_file_to_scratch(
    accession_number: str, group_dir: PosixPath, scratch_dir: PosixPath
):
    if scratch_dir.is_dir() == False:
        scratch_dir.mkdir()

    acc_folder = group_dir / accession_number
    if acc_folder.is_dir():
        acc_number_exists = True
        acc_rna_folder = acc_folder / "RNA" / "FastQ"
        if acc_rna_folder.is_dir():
            rna_folder_exists = True
            for item in acc_rna_folder.iterdir():
                print(item)
                if item.name.endswith(".tar.gz") and item.is_file():

                    fasta_file_exists = True
                    dest_folder = scratch_dir / item.name.split(".tar")[0]
                    if dest_folder.is_dir() == False:
                        dest_folder.mkdir()
                        fastq_folder = dest_folder / "fastq"
                        if fastq_folder.is_dir() == False:
                            fastq_folder.mkdir()
                    dest_file = dest_folder / item.name
                    shutil.copy(str(item), str(dest_file))
                    extract_tar_gz(dest_file, fastq_folder)
                elif item.name.endswith(".fastq.gz") and item.is_file():
                    fasta_file_exists = True
                    dest_folder = scratch_dir / item.name.split(".tar")[0]
                    if dest_folder.is_dir() == False:
                        dest_folder.mkdir()
                        fastq_folder = dest_folder / "fastq"
                        if fastq_folder.is_dir() == False:
                            fastq_folder.mkdir()
                    dest_file = fastq_folder / item.name
                    shutil.copy(str(item), str(dest_file))
            else:
                fasta_file_exists = False
        else:
            rna_folder_exists = False
            fasta_file_exists = False
    else:
        acc_number_exists = False
        rna_folder_exists = False
        fasta_file_exists = False

    return acc_number_exists, rna_folder_exists, fasta_file_exists


def move_rcc_files_to_scratch(
    accession_numbers: list[str],
    group_dir: PosixPath,
    scratch_dir: PosixPath,
):
    acc_numbers_exists = []
    rna_folders_exists = []
    fasta_files_exists = []
    tempus_dirs = []
    for accession_number in accession_numbers:
        for tempus_dir in group_dir.iterdir():
            if tempus_dir.is_dir():
                acc_number_exists, rna_folder_exists, fasta_file_exists = (
                    move_rcc_file_to_scratch(accession_number, tempus_dir, scratch_dir)
                )
                acc_numbers_exists.append(acc_number_exists)
                rna_folders_exists.append(rna_folder_exists)
                fasta_files_exists.append(fasta_file_exists)
                tempus_dirs.append(tempus_dir / accession_number)

    move_report = pd.DataFrame(
        {
            "tempus_dir": tempus_dirs,
            "acc_number_exists": acc_numbers_exists,
            "rna_folder_exists": rna_folders_exists,
            "fasta_file_exists": fasta_files_exists,
        }
    )
    move_report.to_csv(scratch_dir / "move_report.csv", index=False)


if __name__ == "__main__":
    user = sys.argv[1]
    in_dir = sys.argv[2]
    out_dir = sys.argv[3]
    acc_file = pd.read_csv(sys.argv[4])

    move_rcc_files_to_scratch(
        acc_file.accession_id.values,
        Path(f"/group/{user}/{in_dir}"),
        Path(f"/scratch/g/{user}/{out_dir}/"),
    )
