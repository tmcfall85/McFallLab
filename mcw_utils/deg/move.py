import pandas as pd
from pathlib import Path
from datetime import date
import shutil


def main(folder_path, output_folder_path):

    data = {}
    data["accession_id"] = []
    needs_init = True
    genes = []
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
                            dest_file = (
                                Path(output_folder_path)
                                / acc_num
                                / sub_folder.name
                                / fname.name
                            )
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

    args = parser.parse_args()
    main(args.folder_path)
