#!/bin/bash
#SBATCH --job-name=sweeny_rnaseq_test2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --account=tmcfall

module load python/3.9.1
module load rsem/1.3.3
module list
python -m venv rnaseq
source rnaseq/bin/activate
pip install qc-utils==19.8.1 ptools_bin==0.0.7

python rsem_quant.py --anno_bam /home/msochor/mcw_utils/mcw_utils/rna-seq-pipeline/Aligned.toTranscriptome.out.bam --read_strand unstranded --endedness paired --rsem_index /home/msochor/rna-seq-data/ENCFF285DRD.tar.gz --ncpus 12 --ramGB 30
