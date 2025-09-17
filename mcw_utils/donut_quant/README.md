# Running donut quantification experiments

This code will run the donut quantification for McFall Lab experiments.  For details of the methodology please see the published paper (link here: in progress). 

If you want to use this for your project, please be my guest!  I cannot guarantee it'll work out of the box, but hopefully this is a nice starting point for you to do an analysis of this type.  Feel free to email me with questions (msochor@mcw.edu) or make a pull request and ping me.

# Requirements
I am running python 3.12.7.  Please see [requirements](./requirements.txt) for my frozen package versions.  Likely this code will run happily with many versions of packages, but this is a known functional setup.

My GPU is a Nvidia GeForce RTX 4050 with CUDA driver v12.7.33 on a windows laptop.  You could certainly run this without GPU acceleration as its only predict.

## Install mcw_utils
This is needed to import functions within other scripts/jupyter notebooks!

Clone this repo somewhere on your computer, navigate to the top level directory (aka: ./McFallLab) and pip install:

`pip install .`

or editable pip install if you want to play around with the code for your own purposes (this is what I do):

`pip install -e .`

Note: There is a lot of my work all throw into this repo and installed into the python package `mcw_utils` and I use conda to manage different dependencies for different parts of this repo, as different projects need different packages.  So, the above requirements are correct for this folder of the project

# Example
If you want to tldr this README then...

There is an example experiment including expected folder layout (described in detail below) included in this directory.

To run this example, open `run_example.ipynb` in jupyter.  Results are written to the `./example` folder as csv and pkl files.  You better have pip installed `mcw_utils` first!  

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

To run the code from command line using timepoints as normalization:
```python
python donut_quant.py --folder-path base_folder_path
```

To run the code from command line using individual wells as normalization:
```python
python donut_quant.py --folder-path base_folder_path --within-well
```

To run the code from 

It will write csvs for each image timepoint in the basepath as a csv.

# Cropping
This project assumes 96 well plate (12 columns, 8 rows) and because I am lazy, I use PIL to grayscale and crop my images into the actual files that get encoded by ResNet-50.

Your resolution may vary than what I get.  You should use your highest resolution images you can!  

The parameters that are included are:
l: left buffer - essentially the number of pixels to the left of the first well
u: upper buffer - the number of pixels above the top of the first well
w: width - the width of the cropping box, typically the number of pixels wide/tall your well is
step: the number of pixels that needs to be stepped to get to the next well

The default values for this line up with the highest resolution images in the experiments being performed in the McFallLab with our plates.  In the example, b/c I don't want to upload 43MB images, I have significantly smaller images so in the example notebook I use: `l=6, u=6, w=41, step=56`

If you turn on plotting with `show_plot=True` then you can see the cropping.

How you want to crop is up to you, you can crop just the interior of the well if that catches your cells.  Or you can crop the whole well if cells are drifting.  There isn't a perfectly correct answer, its what works for your experiment and setup.

# Normalization - aka within well vs within timepoint
In our experimental setup, the cells are magnetically stabilized into a donut shape (other shapes are possible, and this analysis would still be valid) and then dosed either with a vehicle (no drug) or drug in various doses for a certain amount of time, and then imaged.  Cells are left in an incubator in between, but there is no magnetic field anymore to retain the shape.

Essentially what this means is that the donuts are falling apart even with just vehicle.  What we are attempting to measure is degradation being sped up by drug or different mutational status.

Normalization occurs when we decide what images to compare to as a reference of "not yet degraded".

## Approach 1: comparing within well
This is the most straightforward.  Take timepoint zero (before drug added) for each well as your reference encoding vector.  There will be one reference vector for each well location.

Compare each future timepoint to this reference.

What this tells you is the degradation of each well individually.  Technical replicates can be averaged (and are by default in the code).

### Pros
1. Pretty straightforward to understand
2. Gives you a time series
3. If individual wells are "weird" then you compare starting weird versus it getting weirder over time.  This means its a little more resiliant to things like, I pipetted too hard and disrupted the cells in this particular well.
4. You can use the time series plot to determine where cells are too baseline degraded such that the data is useless (you should also do this qualitatively by looking at images but its better to have both qualitative and quantitative).

## Cons
1. You get no measure at your baseline (because it is the reference)
2. You are measuring the degradation of the cells directly.  So drug effect is the delta of these two and it can be noisy when you take the difference.  You might not see anything

## Approach 2: comparing within timepoints
This approach uses the vehicle controls at each timepoint as the reference vector.  You will have one reference vector for each timepoint.

This tells you how degraded drug is versus vehicle directly.

### Pros
1. You measure drug related degradation directly.  Far less noisy
2. The can use a time series plot to determine where drug is maximally effective.

### Cons
1. You **cannot** plot this as a time series to determine degradation.  Each time should be analyzed individually and treated like an experimental value (i.e. I let the cells have the drug for 12 hours and imaged...).  
2. This approach is susceptible to one weird well throwing off the average.  This can be mitigated by removing that well as a technical replicate (assuming qualitatively you see an experimental error like pipetting disrupting the donut shape)

## Normalization summation
Both approaches are useful.  Run both, the code is pretty quick.  You can answer different questions with different normalization strategies.  

Basically what I'm telling you is to think.  You are a scientist.  Act like it!  Have fun!  Think about what you are doing and what the code is doing.

Please.  Think.  Don't just run code.  This is 2025.  Dark things happen when people just run code.  Or stupid things.