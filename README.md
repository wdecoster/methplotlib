# methylation_browser

This script takes data from nanopolish as processed by calculate_methylation_frequency.py to generate a genome browser view on a window.

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
