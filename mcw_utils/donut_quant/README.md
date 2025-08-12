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