#!/bin/bash
#SBATCH --job-name=sweeny_rnaseq_test
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --account=tmcfall

module load star/2.5.2
module load python/3.9.1
module list
python -m venv rnaseq
source rnaseq/bin/activate
pip install qc-utils==19.8.1 ptools_bin==0.0.7

python align.py --fastqs_R1 /home/msochor/sweeny_data/2_R1_001.fastq.gz --fastqs_R2 /home/msochor/sweeny_data/2_R2_001.fastq.gz --endedness paired --index /home/msochor/rna-seq-data/ENCFF598IDH.tar.gz --bamroot ./outdata/bamtest --ncpus 12 --ramGB 30
