import pytest

import pyranges as pr

import pandas as pd

from methplotlib.differential.differential import main
# def old_method(a_range, b_range, bed):

#     from scipy.stats import fisher_exact

#     def meth_counts(interval):
#         if not interval.empty:
#             all, methylated = interval.as_df()[['calls', "methylated"]].sum()
#             return [all - methylated, methylated]
#         else:
#             return [0, 0]

#     rowdicts = []

#     for chr, begin, end in bed.merge().as_df().itertuples(index=False, name=None):
#         [d1, n1], [d2, n2] = meth_counts(a_range[chr, begin:end]), meth_counts(b_range[chr, begin:end])
#         rowdicts.append({"d1": d1, "d2": d2, "n1": n1, "n2": n2})
#         # print(c1, c2)
#         # print(a, b, c, d)
#         # odr, pv = fisher_exact([c1, c2])
#         # print(odr, pv)
#         # rowdicts.append({"Chromosome": chr, "Start": begin, "End": end, "OR": odr, "P": pv})

#     return pd.DataFrame.from_dict(rowdicts)


columns = ["Chromosome", "Start", "End", "calls", "methylated"]

@pytest.fixture
def pr1():

    df = pd.read_table("tests/d1.tsv.gz", header=0)
    df.columns = columns

    return pr.PyRanges(df)


@pytest.fixture
def pr2():

    df = pd.read_table("tests/d2.tsv.gz", header=0)
    df.columns = columns

    return pr.PyRanges(df)

@pytest.fixture
def bed():

    return pr.read_bed("tests/chr21.bed.gz")

@pytest.fixture
def expected_result():
    df = pd.read_table("tests/result_diff.tsv.gz")
    return df


def test_new_vs_old(pr1, pr2, bed, expected_result):

    result = main(pr1, pr2, bed)
    print(result)
    print(expected_result.head())
    print(expected_result.tail())
    assert 0

