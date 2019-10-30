import pandas as pd
from argparse import ArgumentParser
from statsmodels.sandbox.stats import multicomp


def main():
    args = get_args()
    df = pd.read_csv(args.test_result, sep="\t")
    df["padj"] = bhCorrection(df["p-value"])
    print(df.sort_values(by="padj").to_csv(sep="\t", index=False))


def bhCorrection(s):
    """
    Benjamini-Hochberg correction for a Series of p-values.
    """
    s = s.fillna(1.)
    q = multicomp.multipletests(s, method='fdr_bh')[1][:len(s)]
    q = pd.Series(q[:len(s)], s.index, name='p_adj')
    return q


def get_args():
    parser = ArgumentParser(
        description="Perform multiple testing correction and sort.")
    parser.add_argument("test_result",
                        help="File with results of differential modification "
                             "or allele-specific modification test.")
    return parser.parse_args()


if __name__ == '__main__':
    main()
