# How to run

You need to download 3 files locally, merged data from `tempus_pdac` folder


# to run

```python
from mcw_utils.lyttle_study import label

merged_tempus_fname = 'output_3_19_2025/tempus_json_pdac_rcc_merged.csv'
patients_fname = '/mnt/c/Users/msochor/Downloads/Serum_2023-04-06 - Clinical Data_FINAL.xlsx'
pdac_sotb_fname = '/mnt/c/Users/msochor/Downloads/PDAC_SOTB_2025_03_04_unfiltered.dta'

label(merged_tempus_fname, patients_fname, pdac_sotb_fname)
```