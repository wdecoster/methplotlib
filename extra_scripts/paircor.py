from argparse import ArgumentParser
import pandas as pd


def main():
    args = get_args()
    df = pd.concat([parse_meth_file(f) for f in args.files], axis='columns')


def parse_meth_file(f):
    df = pd.read_csv(f, sep="\t", header=None, names=['chrom', 'pos', 'freq']).drop_duplicates()
    df.index = df['chrom'].astype(str).str.cat(df['pos'].astype(str), sep="_")
    return df.drop(columns=['chrom', 'pos']).rename(columns={'freq': f})


def get_args():
    parser = ArgumentParser(
        description="Make a correlation plot of two results of methylation frequency")
    parser.add_argument("files", help="tsv files with chr, begin, freq (without header)", nargs=2)
    return parser.parse_args()


if __name__ == '__main__':
    main()
