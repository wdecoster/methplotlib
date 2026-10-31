
import pandas as pd
import pyranges as pr
import numpy as np


def fisher_exact(tp, fp, fn, tn, pseudocount=0):
    """
    The same as pyranges.stats.fisher_exact, which passes 64 bit integers
    to the fisher package, while its versions since 0.2.0 require 32 bit integers
    """
    from fisher import pvalue_npy

    tp, fp, fn, tn = (np.asarray(x, dtype=np.uint32) for x in (tp, fp, fn, tn))
    left, right, twosided = pvalue_npy(tp, fp, fn, tn)
    OR = ((tp + pseudocount) / (fp + pseudocount)) / ((fn + pseudocount) / (tn + pseudocount))
    return pd.DataFrame({"OR": OR, "P": twosided, "PLeft": left, "PRight": right})


def _methylated_and_freq_to_zero(df):

    only_in_bed = (df.methylated == -1) & (df.methylated_b == -1)
    df = df[~only_in_bed]

    a_zero = df.methylated == -1
    df.loc[a_zero, ["methylated", "calls"]] = 0
    b_zero = df.methylated_b == -1
    df.loc[b_zero, ["methylated_b", "calls_b"]] = 0

    return df


def count_total_and_methylated(df):
    """For each region in bed, get sum of calls and methylated"""

    grpby = df.groupby("ID")
    return pd.DataFrame({"Chromosome": df.Chromosome.iloc[0],
                         "Start": grpby.Start.first(),
                         "End": grpby.End.first(),
                         "ID": grpby.ID.first(),
                         "calls": grpby.calls.sum(),
                         "methylated": grpby.methylated.sum()})


def merge_regions_with_bed(bed, a, b):

    a = bed.join(a, how="left").apply(count_total_and_methylated).drop(like="_b")
    b = bed.join(b, how="left").apply(count_total_and_methylated).drop(like="_b")

    m = a.join(b, how="outer")

    return m


def main(a, b, bed):
    """1. Find the regions in bed that overlaps either a or/and b.
2. Sum methylations over regions in a/b that overlap bed.
3. Do fisher_exact on the methylation frequencies."""

    bed.ID = np.arange(len(bed))

    if "Strand" in bed:
        bed = bed.drop("Strand")

    m = merge_regions_with_bed(bed, a, b)

    m = m.apply(_methylated_and_freq_to_zero).drop(like="(Start|End|ID)_b|ID")

    m1, c1, m2, c2 = m.methylated, m.calls, m.methylated_b, m.calls_b,
    fe = fisher_exact(c1 - m1, m1, c2 - m2, m2, pseudocount=0.01)
    fe.insert(fe.shape[1], "FDR", pr.stats.fdr(fe.P))

    m = m.insert(fe[['OR', 'P', 'FDR']])

    return m
