[![Build Status](https://travis-ci.com/wdecoster/methplotlib.svg?branch=master)](https://travis-ci.com/wdecoster/methplotlib)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/methplotlib/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/methplotlib/badges/version.svg)](https://anaconda.org/bioconda/methplotlib)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/methplotlib/badges/license.svg)](https://anaconda.org/bioconda/methplotlib)

# METHPLOTLIB

This script generates a browser view on a window using data from  
i) [nanopolish](https://github.com/jts/nanopolish), either as methylation calls or methylation frequencies (as processed by calculate_methylation_frequency.py). The methylation calls can additionally be phased using scripts/annotate_calls_by_phase.and scripts/split_calls_by_phase.py  
ii) [nanocompore](https://github.com/tleonardi/nanocompore)  
iii) in ont-cram format with MM/MP tags according to the SAM specifications and converted using e.g. [this script](https://github.com/kpalin/gcf52ref/blob/f5_to_usam/scripts/extract_methylation_fast5_to_sam.py)  
iv) in bedgraph format  

## INSTALLATION
`pip install methplotlib`

## USAGE
```
methplotlib [-h] [-v] -m METHYLATION [METHYLATION ...] -n NAMES
                   [NAMES ...] -w WINDOW [-g GTF] [-b BED] [-f FASTA]
                   [--simplify] [--split] [--static STATIC] [--smooth SMOOTH]
                   [--dotsize DOTSIZE] [--example] [-o OUTFILE] [-q QCFILE]

plotting nanopolish methylation calls or frequency

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version and exit.
  -m, --methylation METHYLATION [METHYLATION ...]
                        data in nanopolish, nanocompore, ont-cram or bedgraph
                        format
  -n, --names NAMES [NAMES ...]
                        names of datasets in --methylation
  -w, --window WINDOW   window (region) to which the visualisation has to be restricted
  -g, --gtf GTF         add annotation based on a gtf file
  -b, --bed BED         add annotation based on a bed file
  -f, --fasta FASTA     required when --window is an entire chromosome, contig or transcript
  --simplify            simplify annotation track to show genes rather than transcripts
  --split               split, rather than overlay the methylation tracks
  --static              Make a static image of the browser window (filename)
  --binary              Make the nanopolish plot ignorning log likelihood nuances
  --smooth              Rolling window size for averaging frequency values (int)
  --dotsize             Control the size of dots in the per read plots (int)
  --example             Show example command and exit.
  -o, --outfile OUTFILE File to write results to. Default:
                        methylation_browser_{chr}_{start}_{end}.html. Use
                        {region} as a shorthand for {chr}_{start}_{end} in the
                        filename. Missing paths will be created.
  -q, --qcfile QCFILE   File to write the qc report to. Default: The path in
                        outfile prefixed with qc_, default is qc_report_methyl
                        ation_browser_{chr}_{start}_{end}.html. Use {region}
                        as a shorthand for {chr}_{start}_{end} in the
                        filename. Missing paths will be created.

```

## Snakemake workflow
For streamlining nanopolish a Snakefile is included (using snakemake). The workflow uses a config file, of which an example is in this repository.

## Example data
The `examples` folder contains calls and frequencies for the human ACTB gene from PromethION sequencing of NA19240. An example command is available.

## Companion scripts
The `scripts` folder contains scripts for phasing modification calls in haplotypes based on [WhatsHap](https://whatshap.readthedocs.io/en/latest/) phasing, allele specific modification testing for phased data and differential modification testing across subjects.

## TO DO - CONTRIBUTIONS WELCOME
- Outlier detection (in windows) across samples
