import pandas as pd
from argparse import ArgumentParser
import pyranges as pr
from pyranges import PyRanges
from scipy.stats import fisher_exact
import numpy as np


def main():
    args = get_args()
    a_range = pr.concat([pyranges_from_csv(f) for f in args.Agroup])
    b_range = pr.concat([pyranges_from_csv(f) for f in args.Bgroup])

    print("{}\t{}\t{}\t{}\t{}".format("chromosome", "begin", "end", "odds_ratio", "p-value"))
    for chr, begin, end in pr.read_bed(args.bed).merge().as_df().itertuples(index=False, name=None):
        odr, pv = fisher_exact([meth_counts(a_range[chr, begin:end]),
                                meth_counts(b_range[chr, begin:end])])
        if not np.isnan(odr):
            print("{}\t{}\t{}\t{}\t{}".format(chr, begin, end, odr, pv))


def pyranges_from_csv(inputfile):
    colnames = ["Chromosome", "Start", "End", "num_motifs", "calls", "methylated", "freq", "seq"]
    return PyRanges(pd.read_csv(inputfile, sep="\t", names=colnames, header=0))


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
    parser.add_argument("-A", "--Agroup", nargs='+', help="Frequencies of group A.")
    parser.add_argument("-B", "--Bgroup", nargs='+', help="Frequencies of group B.")
    return parser.parse_args()


if __name__ == '__main__':
    main()
