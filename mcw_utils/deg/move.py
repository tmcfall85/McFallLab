import pandas as pd
from pathlib import Path
from datetime import date
import shutil
import os


def main(folder_path, output_folder_path):

    data = {}
    data["accession_id"] = []
    search_dir = Path(folder_path)
    rsem_name = "Aligned.toTranscriptome.out_rsem.genes.results"
    for acc_dir in search_dir.iterdir():
        if acc_dir.is_dir():
            acc_num = acc_dir.stem
            for sub_folder in acc_dir.iterdir():
                if sub_folder.is_dir() and sub_folder.name == "rsem_output":
                    for fname in sub_folder.iterdir():
                        if fname.is_file() and fname.name == rsem_name:
                            print("File to be processed:", fname)
                            dest_path = (
                                Path(output_folder_path) / acc_num / sub_folder.name
                            )
                            dest_path.mkdir(parents=True, exist_ok=True)
                            dest_file = dest_path / fname.name
                            shutil.copy(str(fname), str(dest_file))
                            print("File copied to:", dest_file)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract TPM values for specific genes from RSEM output."
    )
    parser.add_argument(
        "folder_path",
        type=str,
        help="Path to the folder containing subfolders with RSEM output.",
    )
    parser.add_argument(
        "output_folder_path",
        type=str,
        help="Path to the folder where output files will be saved.",
    )

    args = parser.parse_args()
    main(args.folder_path, args.output_folder_path)
