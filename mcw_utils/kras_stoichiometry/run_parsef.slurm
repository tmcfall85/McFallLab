#!/bin/bash
#SBATCH --job-name=measure-ras-parsef
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
if [ -e ./star_output/Aligned.toTranscriptome.out.bam.bai ]; then
    echo "Sorted and indexed BAM file already exists"
else
    echo "Sorting and indexing star alignment file"
    samtools sort ./star_output/Aligned.toTranscriptome.out.bam > ./star_output/Aligned.toTranscriptome.sorted.out.bam
    samtools index ./star_output/Aligned.toTranscriptome.sorted.out.bam
fi
python /home/msochor/McFallLab/mcw_utils/kras_stoichiometry/measure_parsef.py --fname_sorted ./star_output/Aligned.toTranscriptome.sorted.out.bam --fname_seq /scratch/g/tmcfall/rsem.seq --ras kras --out_dir parsef_output
python /home/msochor/McFallLab/mcw_utils/kras_stoichiometry/measure_parsef.py --fname_sorted ./star_output/Aligned.toTranscriptome.sorted.out.bam --fname_seq /scratch/g/tmcfall/rsem.seq --ras hras --out_dir parsef_output