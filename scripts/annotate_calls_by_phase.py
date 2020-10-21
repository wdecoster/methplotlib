import pysam
from argparse import ArgumentParser
import gzip
import sys


def main():
    args = get_args()
    current_chrom = ""
    chrom_seen = []
    meth = gzip.open(args.methylation, 'rt')
    header = next(meth).rstrip().split('\t')
    header.extend(['PS', 'HP'])
    print('\t'.join(header))
    for pos in meth:
        line = pos.rstrip().split('\t')
        if line[0] != current_chrom:
            if line[0] in chrom_seen:
                sys.stderr.write("WARNING: this script is for chromosome-sorted meth files only!")
            sys.stderr.write(f"Switching to {line[0]}\n")
            phased_reads = get_phase_status_dict(bam=args.bam, chrom=line[0])
            current_chrom = line[0]
            chrom_seen.append(line[0])
        line.extend(phased_reads.get(line[4], ["NaN", "NaN"]))
        print('\t'.join(line))


def get_phase_status_dict(bam, chrom):
    return {read.query_name: [str(read.get_tag('PS')), str(read.get_tag('HP'))]
            for read in pysam.AlignmentFile(bam).fetch(contig=chrom)
            if read.has_tag('PS')}


def get_args():
    parser = ArgumentParser(description="Split a nanopolish call-methylation file by haplotypes")
    parser.add_argument("methylation", help="File created by nanopolish call-methylation")
    parser.add_argument("bam", help="bam file created by whatshap haplotag or longshot")
    return parser.parse_args()


if __name__ == '__main__':
    main()
