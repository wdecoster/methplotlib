import pandas as pd
from argparse import ArgumentParser
import pyranges as pr
from pyranges import PyRanges
from scipy.stats import fisher_exact
import numpy as np


def main():
    args = get_args()
    h1 = pyranges_from_csv(args.methylation[0])
    h2 = pyranges_from_csv(args.methylation[1])
    print("{}\t{}\t{}\t{}\t{}".format("chromosome", "begin", "end", "odds_ratio", "p-value"))
    for chr, begin, end in pr.read_bed(args.bed).merge().as_df().itertuples(index=False, name=None):
        odr, pv = fisher_exact([meth_counts(h1[chr, begin:end]),
                                meth_counts(h2[chr, begin:end])])
        if not np.isnan(odr):
            print("{}\t{}\t{}\t{}\t{}".format(chr, begin, end, odr, pv))


def pyranges_from_csv(inputfile):
    colnames = ["Chromosome", "Start", "End", "num_motifs", "calls", "methylated"]
    usecols = [0, 1, 2, 3, 4, 5]
    return PyRanges(pd.read_csv(inputfile, sep="\t", names=colnames, header=0, usecols=usecols))


def meth_counts(interval):
    if not interval.empty:
        all, methylated = interval.as_df()[['calls', "methylated"]].sum()
        return [all - methylated, methylated]
    else:
        return [0, 0]


def get_args():
    parser = ArgumentParser(
        description="Check for modification differences using fisher exact test.")
    parser.add_argument(
        "-b", "--bed", help="Bed file to aggregate modifications on.", required=True)
    parser.add_argument("methylation", nargs=2, help="Haplotype specific frequency files.")
    return parser.parse_args()


if __name__ == '__main__':
    main()
