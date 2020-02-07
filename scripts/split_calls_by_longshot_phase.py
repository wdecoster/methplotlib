import pysam
from argparse import ArgumentParser
import gzip


def main():
    args = get_args()
    h1_reads = set([read.query_name for read in pysam.AlignmentFile(args.bamH1)])
    h2_reads = set([read.query_name for read in pysam.AlignmentFile(args.bamH2)])
    outputs = {f: gzip.open(f"{args.prefix}_{f}_meth.tsv.gz", 'wt') for f in ["H1", "H2", "U"]}
    meths = gzip.open(args.methylation, 'rt')
    header = next(meths)
    for o in outputs.values():
        o.write(header)
    for line in meths:
        outputs[find_read_in_sets(line.split('\t')[4], h1_reads, h2_reads)].write(line)


def find_read_in_sets(string, set_h1, set_h2):
    if string in set_h1:
        return "H1"
    elif string in set_h2:
        return "H2"
    else:
        return "U"


def get_args():
    parser = ArgumentParser(description="Split a nanopolish call-methylation file by haplotypes")
    parser.add_argument("methylation", help="File created by nanopolish call-methylation")
    parser.add_argument("--bamH1", help="bam file created by whatshap haplotag", required=True)
    parser.add_argument("--bamH2", help="bam file created by whatshap haplotag", required=True)
    parser.add_argument("--prefix", help="prefix for the output files", required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()
