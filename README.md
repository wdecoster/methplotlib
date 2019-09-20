# METHPLOTLIB

This script takes data from nanopolish, either as methylation calls or methylation frequencies (as processed by calculate_methylation_frequency.py) to generate a genome browser view on a window.

## INSTALLATION
`pip install methplotlib`

## USAGE
```
methplotlib [-h] -m METHYLATION [METHYLATION ...] -n NAMES [NAMES ...]
                   -w WINDOW [-g GTF] [--simplify] [--split] [--smooth SMOOTH]

Arguments:
  -h, --help            show this help message and exit
  -m METHYLATION [METHYLATION ...], --methylation METHYLATION [METHYLATION ...]
                        nanoplish calls or frequency output
  -n NAMES [NAMES ...], --names NAMES [NAMES ...]
                        names of datasets in --methylation
  -w, --window WINDOW   window (region) to which the visualisation has to be restricted
  -g GTF, --gtf GTF     add annotation based on a gtf file
  --simplify            simplify annotation track to show genes rather than transcripts
  --split               split, rather than overlay the methylation tracks
  --smooth SMOOTH       Smoothen the datapoints, but reduce the details (integer, default=5)
```


## TO DO - CONTRIBUTIONS WELCOME
- Differential methylation analysis (in windows) across groups of samples
- Outlier detection (in windows) across samples
