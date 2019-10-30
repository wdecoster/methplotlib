[![Build Status](https://travis-ci.com/wdecoster/methplotlib.svg?branch=master)](https://travis-ci.com/wdecoster/methplotlib)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/methplotlib/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/methplotlib/badges/version.svg)](https://anaconda.org/bioconda/methplotlib)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/methplotlib/badges/license.svg)](https://anaconda.org/bioconda/methplotlib)

# METHPLOTLIB

This script takes data from nanopolish, either as methylation calls or methylation frequencies (as processed by calculate_methylation_frequency.py) to generate a genome browser view on a window.

## INSTALLATION
`pip install methplotlib`

## USAGE
```
methplotlib [-h] [-v] -m METHYLATION [METHYLATION ...] -n NAMES
               [NAMES ...] -w WINDOW [-g GTF] [-b BED] [--simplify]
               [--split] [--smooth SMOOTH]

Arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version and exit.
  -m, --methylation METHYLATION [METHYLATION ...]
                        nanoplish calls or frequency output
  -n, --names NAMES [NAMES ...]
                        names of datasets in --methylation
  -w, --window WINDOW   window (region) to which the visualisation has to be restricted
  -g, --gtf GTF         add annotation based on a gtf file
  -b, --bed BED         add annotation based on a bed file matching to your reference genome
  --simplify            simplify annotation track to show genes rather than transcripts
  --split               split, rather than overlay the methylation frequency tracks
  --smooth SMOOTH       Smoothen the datapoints of frequencies, but reduce the details (integer, default=5)
```

## Snakemake workflow
For streamlining nanopolish a Snakefile is included (using snakemake). The workflow uses a config file, of which an example is in this repository.

## Test data
The `test` folder contains calls and frequencies for the human ACTB gene from PromethION sequencing of NA19240. An example command is available.

## Companion scripts
The `scripts` folder contains scripts for phasing modification calls in haplotypes based on [WhatsHap](https://whatshap.readthedocs.io/en/latest/) phasing, allele specific modification testing for phased data and differential modification testing across subjects.

## TO DO - CONTRIBUTIONS WELCOME
- Outlier detection (in windows) across samples
