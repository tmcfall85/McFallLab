from pathlib import Path
from ffile import Ffile
import sys


def pipeline(user, out_dir):
    endedness = "single"
    scratch_dir = Path(f"/scratch/g/{user}/{out_dir}/")
    star_version = "2.5.2"
    rsem_version = "1.3.3"
    python_version = "3.9.1"
    git_branch = "issue_1"
    cpus_per_task = 12
    ram_gb_per_task = 30
    star_index = "/home/msochor/rna-seq-data/ENCFF598IDH.tar.gz"
    rsem_index = "/home/msochor/rna-seq-data/ENCFF285DRD.tar.gz"
    star_template = Ffile("run_align.slurm.template")
    rsem_template = Ffile("run_rsem_quant.slurm.template")
    for subfolder in scratch_dir.iterdir():
        if subfolder.is_dir():
            fastqs = []
            for item in subfolder.iterdir():
                if item.suffix == ".gz":
                    fastqs.append(item)
            if len(fastqs) == 2:
                fastqs_r1 = None
                fastqs_r2 = None
                for fastq in fastqs:
                    if fastq.lower().find("r1") > -1:
                        fastqs_r1 = fastq
                    elif fastq.lower().find("r1") > -1:
                        fastqs_r2 = fastq

                if fastqs_r1 == None or fastqs_r2 == None:
                    print("cannot assign fastqs to r1 and r2!!!")
                    print(fastqs)
                    print("doing it randomly, I guess")
                    fastqs_r1 = fastq[0]
                    fastqs_r2 = fastq[1]
            elif len(fastqs) == 1:
                fastqs_r1 = fastqs[0]
                fastqs_r2 = "none"
            else:
                print("Too many or too few fastqs found!!")
                print(fastqs)
                raise NameError
            fastqs_r1 = fastqs[0]
            with open(subfolder / "run.slurm", "w") as fp:
                fp.write(
                    star_template.f(
                        fname=item.stem.split(".")[0],
                        cpus_per_task=cpus_per_task,
                        account=user,
                        star_version=star_version,
                        python_version=python_version,
                        rsem_version=rsem_version,
                        git_branch=git_branch,
                        fastqs_r1=fastqs_r1,
                        fastqs_r2=fastqs_r2,
                        endedness=endedness,
                        star_index=star_index,
                        rsem_index=rsem_index,
                        ram_gb_per_task=ram_gb_per_task,
                    )
                )
            """
            with open(subfolder / "run_rsem_quant.slurm", "w") as fp:
                fp.write(
                    rsem_template.f(
                        fname=item.stem.split(".")[0],
                        cpus_per_task=cpus_per_task,
                        account=user,
                        python_version=python_version,
                        git_branch=git_branch,
                        anno_bam=anno_bam,
                        endedness=endedness,
                        ram_gb_per_task=ram_gb_per_task,
                    )
                )
            """


if __name__ == "__main__":
    user = sys.argv[1]
    out_dir = sys.argv[2]

    pipeline(user, out_dir)
