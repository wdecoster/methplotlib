from argparse import ArgumentParser
import gzip


def main():
    args = get_args()
    meth = gzip.open(args.input, 'rt')
    next(meth)  # skip the header
    with gzip.open(args.output, 'wt') as out:
        for line in meth:
            ll = line.split('\t')
            out.write(f'{ll[0]}:{ll[1]}-{ll[2]}\t{convert_OR(ll[7])}\n')
    meth.close()


def convert_OR(odr):
    return odr if float(odr) > 1 else str(1/float(odr))


def get_args():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input",
                        help="gzipped input file from allele_specific_methylation",
                        required=True)
    parser.add_argument("-o", "--output", help="gzipped output file", required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()
