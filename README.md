# METHPLOTLIB

This script takes data from nanopolish as processed by calculate_methylation_frequency.py to generate a genome browser view on a window.

## INSTALLATION
`pip install methplotlib`

## USAGE
```
methylation_frequency_browser.py [-h] -m METHYLATION [METHYLATION ...]
                                        -n NAMES [NAMES ...] -w WINDOW
                                        [-g GTF]


optional arguments:
  -h, --help            show this help message and exit
  -m, --methylation METHYLATION [METHYLATION ...]
                        output of calculate_methylation_frequency.py
  -n, --names NAMES     names of datasets in --methylation
  -w, --window WINDOW   window (region) to which the visualisation has to be restricted e.g. chr7:12345-23456
  -g, --gtf GTF     add annotation based on a gtf file matching to your
                        reference genome
```


## TO DO - CONTRIBUTIONS WELCOME
- Make Region() accept coordinate notation with comma's such as  chr14:105,476,219-105,479,218
- Adapt to also use methylation data from tombo
- Add plot showing (pairwise) correlation between methylation results
- Play around with parameter of sliding window average and find an optimal default
- Add snakemake workflow to streamline preprocessing of nanopolish: from fast5 directory to methylation frequency
- Add a mode to plot raw methylation signals rather than summarized frequencies
- Think about how to solve (split?) a region/window which is too large to plot in one file
- Speed up generating plots of annotation
- Simplify annotation plotting - avoiding transcripts which are duplicates
- Differential methylation analysis (in windows) across groups of samples
- Outlier detection (in windows) across samples
