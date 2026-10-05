Run differential expression on PDAC sotb data living in both RCC and outside of RCC:

# RCC stuff
Current RSEM results are stored on RCC: `/group/dseo/work/tempus_rcc_results`

0. Log into RCC, open a terminal and navigate to this folder.  Load up python with:
   
`module load python/3.9.1`

1. Copy files to scratch

```python
python move.py /group/dseo/work/tempus_rcc_results /scratch/g/tmcfall/tempus_rcc_results
```

2. Concatenate TPM results from RSEM for each gene

```python
python quantify_genes.py /scratch/g/tmcfall/tempus_rcc_results
```

cleanup files in scratch please and thank you

`rm -rf /scratch/g/tmcfall/tempus_rcc_results`

3. Download results to your computer, should be `deg_gene_expression_results_DATE.csv`
   
# Local Stuff

Locally you will now run the differential expression analysis.

## Basic mode
This will run DeSeq locally

Usage: python deg.py <deg_file> <dta_file> <merge_file>

## End-to-end K-fold modeling
This mode will split the patients into stratified K-folds and run DeSeq on each training set as well as do the accounting for end-to-end k-fold modeling downstream

Usage: python deg.py <deg_file> <dta_file> <merge_file> <k_fold> <holdout_fraction>