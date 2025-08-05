#!/bin/bash
#SBATCH --job-name=create_pdac_deg
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --account=tmcfall

module load python/3.9.1

python quantify_genes.py /scratch/g/tmcfall/tempus_rsem_results