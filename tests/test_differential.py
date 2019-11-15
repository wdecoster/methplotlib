from io import StringIO

import pytest

import pyranges as pr

import pandas as pd


from methplotlib.differential.differential import main, merge_regions_with_bed

def old_method(a_range, b_range, bed):

    from scipy.stats import fisher_exact

    def meth_counts(interval):
        if not interval.empty:
            all, methylated = interval.as_df()[['calls', "methylated"]].sum()
            return [all - methylated, methylated]
        else:
            return [0, 0]

    rowdicts = []

    for chr, begin, end in bed.merge().as_df().itertuples(index=False, name=None):
        [d1, n1], [d2, n2] = meth_counts(a_range[chr, begin:end]), meth_counts(b_range[chr, begin:end])
        rowdicts.append({"Chromosome": chr, "Start": begin, "End": end, "d1": d1, "d2": d2, "n1": n1, "n2": n2})
        # print(c1, c2)
        # print(a, b, c, d)
        # odr, pv = fisher_exact([c1, c2])
        # print(odr, pv)
        # rowdicts.append({"Chromosome": chr, "Start": begin, "End": end, "OR": odr, "P": pv})

    return pd.DataFrame.from_dict(rowdicts)


columns = ["Chromosome", "Start", "End", "calls", "methylated"]

@pytest.fixture
def pr1():

    df = pd.read_csv("tests/d1.tsv.gz", header=0, sep="\t")
    df.columns = columns

    return pr.PyRanges(df)


@pytest.fixture
def pr2():

    df = pd.read_csv("tests/d2.tsv.gz", header=0, sep="\t")
    df.columns = columns

    return pr.PyRanges(df)

@pytest.fixture
def bed():

    return pr.read_bed("tests/chr21.bed.gz")

@pytest.fixture
def expected_result():
    df = pd.read_csv("tests/result_diff.tsv.gz", sep="\t")
    return df


def test_new_vs_old(pr1, pr2, bed, expected_result):

    result = main(pr1, pr2, bed)
    # old_result = old_method(pr1, pr2, bed)
    counts = result.sort().df[["calls", "methylated", "calls_b", "methylated_b"]]
    counts.loc[:, "calls"] = counts.calls - counts.methylated
    counts.loc[:, "calls_b"] = counts.calls_b - counts.methylated_b
    # old_result.to_csv("expected.tsv", sep="\t")
    # print(old_result)
    assert (counts.calls == expected_result.d1.values).all()
    assert (counts.calls_b == expected_result.d2.values).all()
    assert (counts.methylated == expected_result.n1.values).all()
    assert (counts.methylated_b == expected_result.n2.values).all()


@pytest.fixture()
def simple_bed():

    gr = pr.PyRanges(chromosomes=[1, 1, 1, 1, 1], starts=[0, 10, 20, 30, 100], ends=[5, 15, 25, 35, 105])
    gr.ID = list(range(len(gr)))
    return gr


@pytest.fixture()
def simple_m1():

    m1 = pr.PyRanges(chromosomes=[1, 1, 1, 1], starts=[0, 20, 30, 50], ends=[5, 25, 35, 55])
    m1.calls = [1, 2, 3, 1]
    m1.methylated = [0, 2, 0, 1]
    return m1


@pytest.fixture()
def simple_m2():

    m2 = pr.PyRanges(chromosomes=[1, 1, 1, 1, 1], starts=[0, 10, 30, 40, 55], ends=[5, 15, 35, 45, 55])
    m2.calls = [1, 1, 2, 3, 5]
    m2.methylated = [0, 0, 2, 0, 5]

    return m2

def test_simple(simple_bed, simple_m1, simple_m2):

    print(simple_bed, simple_m1, simple_m2)
    result = merge_regions_with_bed(simple_bed, simple_m1, simple_m2)
    print(result)
    # means that the joins were incorrect 
    assert (result.Start != -1).all()




