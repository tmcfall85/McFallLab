#!/bin/bash
#SBATCH --job-name=measure-ras-parsef
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#SBATCH --account=tmcfall

module load miniforge
module load samtools/1.2
module list

conda activate mcw_utils
mkdir -p /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/parsef
samtools sort /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/star_output/Aligned.to.Transcriptome.out.bam > /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/star_output/Aligned.to.Transcriptome.sorted.out.bam
samtools index /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/star_output/Aligned.to.Transcriptome.sorted.out.bam
python measure_parsef.py --fname_sorted /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/star_output/Aligned.to.Transcriptome.sorted.out.bam --fname_seq /group/tmcfall/work/rsem.seq --ras kras --out_dir /scratch/g/dseo/tempus_rna_seq_part3/TL-19-043DB3/parsef