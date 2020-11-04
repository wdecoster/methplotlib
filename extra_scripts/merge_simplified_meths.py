from argparse import ArgumentParser
import pandas as pd


def main():
    args = get_args()
    print(pd.concat([pd.read_csv(f,
                                 sep="\t",
                                 index_col=0,
                                 header=None,
                                 names=['locus', f.replace('_haplotype_specific_meth.tsv.gz', '')])
                     for f in args.files],
                    axis='columns')
          .to_csv(sep="\t", index_label='locus', na_rep='NaN')
          )


def get_args():
    parser = ArgumentParser()
    parser.add_argument("files", nargs="+", help="simplified haplotype specific meth files")
    return parser.parse_args()


if __name__ == '__main__':
    main()
