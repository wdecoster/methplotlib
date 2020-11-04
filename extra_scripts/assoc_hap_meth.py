from argparse import ArgumentParser
import pandas as pd
import re
import sys
from fisher import pvalue


def main():
    args = get_args()
    df = pd.read_csv(args.meth, sep="\t", index_col=0)
    df.columns = [re.sub('_v4.*', '', d) for d in df.columns]
    sample_info = pd.read_csv(args.sampleinfo, sep="\t", index_col=0)
    cases, controls = parse_sample_info(sample_info)
    df = df[df.count(axis=1) != 1]
    df = (df > 1).astype(int)
    df_pat = df[cases]
    df_con = df[controls]
    for (locus, *pats), cons in zip(df_pat.itertuples(index=True, name=None),
                                    df_con.itertuples(index=False, name=None)):
        p = pvalue(pats.count(1), pats.count(0), cons.count(1), cons.count(0))
        print(f"{locus}\t{p.two_tail}")


def parse_sample_info(sample_info):
    cases = []
    controls = []
    for sample, group in sample_info.itertuples(index=True, name=None):
        if group in ['FTLD-FUS', 'FTLD-TDP-A', 'FTLD-TDP-C', 'FTLD-TDP-B']:
            cases.append(sample)
        elif group == 'Control':
            controls.append(sample)
        else:
            sys.exit(f"Unexpected group name: {group}")
    return cases, controls


def get_args():
    parser = ArgumentParser()
    parser.add_argument("meth", help="merged haplotype specific methylation file")
    parser.add_argument("sampleinfo", help="group information for included samples")
    return parser.parse_args()


if __name__ == '__main__':
    main()
