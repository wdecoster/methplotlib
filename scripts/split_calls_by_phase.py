from argparse import ArgumentParser
import gzip


def main():
    args = get_args()
    phased_calls = gzip.open(args.phased_methylation, 'rt')

    phase1 = gzip.open(args.prefix + "_phase1.tsv.gz", 'wt')
    phase2 = gzip.open(args.prefix + "_phase2.tsv.gz", 'wt')
    uphase = gzip.open(args.prefix + "_unphased.tsv.gz", 'wt')

    header = next(phased_calls)
    for output in [phase1, phase2, uphase]:
        output.write(header)

    for line in phased_calls:
        if line.endswith('1.0\n'):
            phase1.write(line)
        elif line.endswith('2.0\n'):
            phase2.write(line)
        else:
            uphase.write(line)


def get_args():
    parser = ArgumentParser(description="Split file with phased calls by phase.")
    parser.add_argument("phased_methylation", help="File created by annotate_calls_by_phase.py")
    parser.add_argument("prefix", help="Prefix for output files")
    return parser.parse_args()


if __name__ == '__main__':
    main()
