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
- Adapt to also use methylation data from tombo
- Add plot showing (pairwise) correlation between methylation results
- Add a mode to plot raw methylation signals rather than summarized frequencies
- Differential methylation analysis (in windows) across groups of samples
- Outlier detection (in windows) across samples
- Add pipeline and visualization for phased methylation
- Think about adding a coverage track - maybe based on the called_sites column
