# What this
This will quantify hras, nras, and kras transcripts from a bam file.  For kras, it looks for the G12D mutant and pulls those out separately from the kras count.  i.e. if there are 100 kras transcripts, 12 cover the G12 region, and 8 are WT it should output, 90 kras, 8 WT, 4 G12D as counts

Todo: quantify other mutants

# Getting started
Run the tempus_pdac merge job, and upload whatever merged file has RNA files that you want to quantify (e.g. RNA_rcc_bam_with_pdac_sotb_mrn.csv)

make a scratch dir: (e.g. /scratch/dseo/rna_seq_temp)

update your csv to home dir, them move files as seen below

# Moving files on RCC

From this directory:

update as needed

`python move.py dseo work/tempus rna_seq_temp ~/RNA_rcc_bam_with_pdac_sotb_mrn.csv`

# Run files after they have been moved

update run.slurm if needed, then run:

`sbatch run.slurm`