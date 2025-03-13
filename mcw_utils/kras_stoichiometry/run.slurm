#!/bin/bash
#SBATCH --job-name=measure-ras-freq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --account=dseo

module load python/3.9.1
module load samtools/1.2
module list

python measure.py /scratch/g/dseo/rna_seq_temp