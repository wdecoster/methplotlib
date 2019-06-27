import pysam
from argparse import ArgumentParser
import pandas as pd


def main():
    args = get_args()
    calls = pd.read_csv(args.methylation, sep="\t")
    calls['PS'] = None
    calls['HP'] = None
    for read in pysam.AlignmentFile("haplotagged.bam").fetch():
        if read.has_tag('PS'):
            calls.loc[calls['read_name'] == read.query_name, ["PS", "HP"]] = (
                read.get_tag('PS'), read.get_tag('HP'))
    if args.output == 'stdout':
        print(calls.to_csv(path_or_buf=None, sep="\t", na_rep="NaN", index=False))
    else:
        calls.to_csv(path_or_buf=args.output, sep="\t", na_rep="NaN", index=False)


def get_args():
    parser = ArgumentParser(description="Split a nanopolish call-methylation file by haplotypes")
    parser.add_argument("methylation", help="File created by nanopolish call-methylation")
    parser.add_argument("bam", help="bam file created by whatshap haplotag")
    parser.add_argument("-o", "--output",
                        help="Output file to write to. Default is to stdout.",
                        default='stdout')
    return parser.parse_args()


if __name__ == '__main__':
    main()
