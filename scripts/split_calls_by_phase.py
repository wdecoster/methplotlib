import pysam
from argparse import ArgumentParser
import pandas as pd


def main():
    args = get_args()
    phase_info = pd.DataFrame(
        data=[(read.query_name, read.get_tag('PS'), read.get_tag('HP'))
              for read in pysam.AlignmentFile(args.bam).fetch() if read.has_tag('PS')],
        columns=['read_name', 'PS', 'HP']) \
        .set_index('read_name')
    print(pd.read_csv(args.methylation, sep="\t")
          .join(phase_info, on='read_name')
          .to_csv(path_or_buf=None, sep="\t", na_rep="NaN", index=False))


def get_args():
    parser = ArgumentParser(description="Split a nanopolish call-methylation file by haplotypes")
    parser.add_argument("methylation", help="File created by nanopolish call-methylation")
    parser.add_argument("bam", help="bam file created by whatshap haplotag")
    return parser.parse_args()


if __name__ == '__main__':
    main()
