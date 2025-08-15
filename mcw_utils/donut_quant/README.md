# Running donut quantification experiments

These are a bit bespoke at the start, but hopefully we eventually get to a point where its more generalized

# 6_19_25 experiment
I hand cropped images Raven sent me in GIMP and saved them locally.  

The code can optionally use replicates from vehicle, RMC, or TRM wells.  To include replicates, you basically provide it: `[[1,2],[2],[1,2]]` for example uses vehicle replicates 1 and 2, RMC replicate 2, and TRM replicates 1 and 2

Results saved locally to CSVs

```python
python donut_quant_6_19_25.py '/mnt/c/Users/msochor/Downloads/6.19.25_before drug/0005587_01' [[1,2],[1,2],[1,2]]
python donut_quant_6_19_25.py '/mnt/c/Users/msochor/Downloads/6.19.25_before drug/0005587_01' [[2],[1],[1,2]]
python donut_quant_6_19_25.py '/mnt/c/Users/msochor/Downloads/6.19.25_before drug/0005587_01' [[1],[2],[1,2]]
```

# Previous PDAC experiment
I don't know exactly when Raven ran this one, but she sent me another old experiment.  The resolution was all over the place, and I had to hand crap the major images again.  So this is not standardized, but its closer to it.  

```python
python donut_quant_previous_pdac.py /mnt/c/Users/msochor/Downloads trm_results.csv
```

# General PDAC experiment
This is for a 96 well plate experiment.  It assumes there is a folder:

```
base_path/image_folder1/******/******.tif
base_path/image_folder2/******/******.tif
...
base_path/image_folderN/******/******.tif
base_path/samples.csv
base_path/drug.csv
base_path/files.csv
```

You need to create the csv:
example samples.csv
```
None	None	None	None	None	None	None	None	None	None	None	None
None	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	None
None	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	None
None	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	G12D	None
None	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	None
None	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	None
None	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	G12R	None
None	None	None	None	None	None	None	None	None	None	None	None
```

example drug.csv
```
None	None	None	None	None	None	None	None	None	None	None	None
None	0	10	18	32	56	100	178	316	562	1000	None
None	0	10	18	32	56	100	178	316	562	1000	None
None	0	10	18	32	56	100	178	316	562	1000	None
None	0	10	18	32	56	100	178	316	562	1000	None
None	0	10	18	32	56	100	178	316	562	1000	None
None	0	10	18	32	56	100	178	316	562	1000	None
None	None	None	None	None	None	None	None	None	None	None	None

```

and files.csv
```
time	folder
-1	7.23.25_before_drug
0	7.23.25_after_drug
1	7.23.25_1_hr
12	7.23.25_12_hr
24	7.24.25_24_hr
48	7.25.25_48hrs
72	7.26.25_72_hr
```

where the folder column lets you know the `image_folder1` above.

To run the code:
```python
python donut_quant.py base_folder_path
```

It will write csvs for each image timepoint in the basepath as a csv.