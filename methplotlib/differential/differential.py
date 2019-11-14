
import pandas as pd
import pyranges as pr
from pyranges import PyRanges
import numpy as np


def _methylated_and_freq_to_zero(df):

    a_zero = df.Start == -1
    df.loc[a_zero, ["methylated", "calls"]] = 0
    b_zero = df.Start_b == -1
    df.loc[b_zero, ["methylated_b", "calls_b"]] = 0

    return df


def count_total_and_methylated(df):
    """For each region in bed, get sum of calls and methylated"""

    grpby = df.groupby("ID")
    total = grpby.calls.sum()
    methylated = grpby.methylated.sum()

    chromosome = df.Chromosome.iloc[0]
    start = grpby.Start.first()
    end = grpby.End.first()
    region_id = grpby.ID.first()

    return pd.DataFrame({"Chromosome": chromosome, "Start": start, "End": end, "ID": region_id,
                         "calls": total, "methylated": methylated})


def main(a, b, bed):

    bed.ID = np.arange(len(bed))

    if "Strand" in bed:
        bed = bed.drop("Strand")

    a = bed.join(a).apply(count_total_and_methylated).drop(like="_b")
    b = bed.join(b).apply(count_total_and_methylated).drop(like="_b")

    m = a.join(b, how="outer")

    m = m.apply(_methylated_and_freq_to_zero).drop(like="(Start|End|ID)_b|ID")

    m1, c1, m2, c2 = m.methylated, m.calls, m.methylated_b, m.calls_b,
    fe = pr.stats.fisher_exact(c1 - m1, m1, c2 - m2, m2, pseudocount=0.01, alternative="twosided")

    fe.insert(fe.shape[1], "ORFDR", pr.stats.fdr(fe.P))

    m = m + fe

    return m
