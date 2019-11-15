import pandas as pd
from pyranges import PyRanges


def methylation_pyranges_from_csv(inputfile):
    colnames = ["Chromosome", "Start", "End", "calls", "methylated"]
    return PyRanges(pd.read_csv(inputfile,
                                sep="\t",
                                names=colnames,
                                header=0,
                                usecols=[0, 1, 2, 4, 5]))
