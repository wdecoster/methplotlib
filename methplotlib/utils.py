from argparse import ArgumentParser
import sys


class Region(object):
    def __init__(self, region):
        self.chromosome, interval = region.replace(',', '').split(':')
        self.begin, self.end = [int(i) for i in interval.split('-')]


def get_args():
    parser = ArgumentParser(description="plotting methylation frequency")
    parser.add_argument("-m", "--methylation",
                        nargs='+',
                        help="output of calculate_methylation_frequency.py",
                        required=True)
    parser.add_argument("-n", "--names",
                        nargs='+',
                        help="names of datasets in --methylation",
                        required=True)
    parser.add_argument("-w", "--window",
                        help="window (region) to which the visualisation has to be restricted",
                        required=True)
    parser.add_argument("-g", "--gtf",
                        help="add annotation based on a gtf file matching to your reference genome")
    parser.add_argument("-s", "--smooth",
                        help="Smoothen the datapoints, but reduce the details",
                        type=int,
                        default=5)
    args = parser.parse_args()
    if not len(args.names) == len(args.methylation):
        sys.exit("INPUT ERROR: Expecting the same number of names as datasets!")
    return args
