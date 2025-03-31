from pathlib import Path
from shutil import copy2
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


def main(user: str, subfolder: str, outdir: Path, copy_bam=bool):
    scratch_dir = Path(f"/scratch/g/{user}/{subfolder}")
    if scratch_dir.is_dir() == False:
        print(f"Not a valid directory! {scratch_dir}")
        print("Exiting.")
        return
    if outdir.is_dir() == False:
        outdir.mkdir()
    for tempus_folder in scratch_dir.iterdir():
        if tempus_folder.is_dir():
            if verify_files(tempus_folder):
                print(f"Folder {tempus_folder} has completed data")
                out_tempus = outdir / tempus_folder.stem
                if out_tempus.is_dir():
                    print("Folder is already copied.  Skipping.")
                    print(
                        f"If you think this is not correct, please go check for what is in {out_tempus}"
                    )
                else:
                    out_tempus.mkdir()
                    for file in tempus_folder.iterdir():
                        if (
                            file.is_file()
                            and file.suffix == ".out"
                            and file.name.startswith(tempus_folder.stem)
                        ):
                            print(f"Copying rna seq pipeline log: {file}")
                            copy2(file, outdir / "mwc_utils_rna_seq_pipeline.log")
                        elif file.is_file() and file.name == "run.slurm":
                            print(f"Copying run.slurm: {file}")
                            copy2(file, outdir / "run.slurm")
                    star_output = tempus_folder / "star_output"
                    star_output_dest = out_tempus / "star_output"
                    if star_output.is_dir():
                        star_output_dest.mkdir()
                        print("Copying STAR output files:")
                        for file in star_output.iterdir():
                            if file.suffix == ".bam":
                                if copy_bam:
                                    print(file)
                                    copy2(file, star_output_dest / file.name)
                                else:
                                    print(f"skipping bam file b/c its too swol {file}")
                            else:
                                print(file)
                                copy2(file, star_output_dest / file.name)
                    else:
                        print("MISSING STAR OUTPUT! RUH ROH RAGGY")
                    rsem_output = tempus_folder / "rsem_output"
                    rsem_output_dest = out_tempus / "rsem_output"
                    if rsem_output.is_dir():
                        rsem_output_dest.mkdir()
                        print("Copying rsem output files:")
                        for file in rsem_output.iterdir():
                            print(file)
                            copy2(file, rsem_output_dest / file.name)
                    else:
                        print("MISSING RSEM OUTPUT! RUH ROH RAGGY")

            else:
                print(f"Folder {tempus_folder} is missing data!")


if __name__ == "__main__":
    user = sys.argv[1]
    subfolder = sys.argv[2]
    outdir = Path(sys.argv[3])
    if len(sys.argv) > 4:
        if sys.argv[4].lower() == "--save-bam":
            copy_bam = True
        else:
            copy_bam = False
    else:
        copy_bam = False

    main(user, subfolder, outdir, copy_bam)
