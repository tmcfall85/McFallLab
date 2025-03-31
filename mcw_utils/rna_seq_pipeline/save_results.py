from pathlib import Path
import sys


def verify_files(tempus_folder):
    # basically these files should exist if the process completed successfully
    rsem_output = tempus_folder / "rsem_output"
    if rsem_output.is_dir():
        num_genes = (
            rsem_output / "Aligned.toTranscriptome.out_number_of_genes_detected.json"
        )
        if num_genes.is_file():
            return True

    return False


def main(user: str, subfolder: str):
    scratch_dir = Path(f"/scratch/g/{user}/{subfolder}")
    if scratch_dir.is_dir() == False:
        print(f"Not a valid directory! {scratch_dir}")
        print("Exiting.")
        return
    for tempus_folder in scratch_dir.iterdir():
        if tempus_folder.is_dir():
            if verify_files(tempus_folder):
                print(f"Folder {tempus_folder} has completed data")
            else:
                print(f"Folder {tempus_folder} is missing data!")


if __name__ == "__main__":
    user = sys.argv[1]
    subfolder = sys.argv[2]
    main(user, subfolder)
