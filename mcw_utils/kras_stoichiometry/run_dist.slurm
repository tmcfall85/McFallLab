#!/bin/bash
#SBATCH --job-name=measure-ras-dist
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=31gb
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --account=tmcfall
#SBATCH --partition=bigmem

module load miniforge
module load samtools/1.20
module list

conda activate mcw_utils
python /home/msochor/McFallLab/mcw_utils/kras_stoichiometry/measure_nras_count.py --fname_sorted ./star_output/Aligned.toTranscriptome.sorted.out.bam --fname_seq /scratch/g/tmcfall/rsem.seq --out_dir dist_output
