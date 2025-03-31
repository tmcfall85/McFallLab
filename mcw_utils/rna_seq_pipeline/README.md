# What this
This will run STAR and RSEM on raw fasta files.  It generates BAM alignment files from STAR as well as gene and isoform count files from RSEM.

# Getting started
Basically the process is to move files to /scratch dir on RCC, generate slurm run files for each folder, and then manually kick off each slurm file (this could be automated, but in the interest of not bombarding a shared resource like RCC, I decided to just leave it manual for now).

It typically takes 1-2 hours to do the alignment and feature counting.  RCC can run these in parallel as resources allow.


# Moving files on RCC 

Note: this works for the tempus files stored in the dseo directory.  If there is a different directory structure and/or file structure, it will need to be amended.  The tempus files are 2 fasta files gzipped together, so this copies that file and unzips it on scratch

From this directory on RCC run this:

`python move.py dseo work/tempus tempus_rna_seq ~/RNA_rcc_bam_with_pdac_sotb_mrn.csv`

# Generate slurm files

Slurm files are hardcoded for local files.  They also define what version of each program is run (star, rsem, python, etc)

`python pipeline.py dseo tempus_rna_seq`

# Run each manually

Navigate over to the scratch subfolder and do:

`sbatch run.slurm`

# Copy results back and cleanup scratch

Scratch isn't long term storage, so we need to copy the resultant files back:

`python save_results.py dseo tempus_rna_seq`

Then cleanup the dir we made in scratch.  NOTE: MAKE SURE YOU HAVE EVERYTHING YOU WANT SAVED BEFORE YOU DO THIS OR YOU GOTS TO RUN LOTS OF THINGS AGAIN AND THIS TAKES 