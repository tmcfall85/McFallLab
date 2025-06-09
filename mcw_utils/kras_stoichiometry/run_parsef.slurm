#!/bin/bash
#SBATCH --job-name=measure-ras-parsef
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --account=tmcfall

module load python/3.9.1
module load samtools/1.2
module list

source stoichiometryEnv/bin/activate

python -m venv parsef
source parsef/bin/activate
rm -rf McFallLab/
git clone https://github.com/tmcfall85/McFallLab.git
cd McFallLab 
git checkout issue_11
git rev-parse HEAD
pip install -e .
cd ../
mkdir -p /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/parsef
samtools sort /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/star_output/Aligned.to.Transcriptome.out.bam > /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/star_output/Aligned.to.Transcriptome.sorted.out.bam
samtools index /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/star_output/Aligned.to.Transcriptome.sorted.out.bam
python measure_parsef.py /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/star_output/Aligned.to.Transcriptome.sorted.out.bam /group/tmcfall/work/rsem.seq kras /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/parsef