#!/bin/bash
#SBATCH --job-name=measure-ras-freq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --partition=bigmem
#SBATCH --account=dseo

module load python/3.9.1
module load samtools/1.15.1
module list

python -m venv stoichiometryEnv
source stoichiometryEnv/bin/activate

pip install pysam
pip install pandas

python measure.py /scratch/g/tmcfall/measure_tempus