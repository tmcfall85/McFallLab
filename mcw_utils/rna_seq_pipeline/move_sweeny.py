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


def move_rcc_files_to_scratch(
    group_dir,
    scratch_dir,
):
    if scratch_dir.is_dir() == False:
        scratch_dir.mkdir()
    for item in group_dir.iterdir():

        if item.name.endswith(".fastq.gz") and item.is_file():
            dest_folder = scratch_dir / item.name.split(".tar")[0]
            if dest_folder.is_dir() == False:
                dest_folder.mkdir()
            fastq_folder = dest_folder / "fastq"
            if fastq_folder.is_dir() == False:
                fastq_folder.mkdir()
            dest_file = fastq_folder / item.name
            print(f"Copying {item} to {dest_file}")
            shutil.copy(str(item), str(dest_file))


if __name__ == "__main__":
    user = sys.argv[1]
    in_dir = sys.argv[2]
    out_dir = sys.argv[3]

    move_rcc_files_to_scratch(
        Path(f"/group/{user}/{in_dir}"),
        Path(f"/scratch/g/{user}/{out_dir}/"),
    )
