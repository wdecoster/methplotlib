from argparse import ArgumentParser
import sys
from math import ceil


class Region(object):
    def __init__(self, region):
        self.chromosome, interval = region.replace(',', '').split(':')
        self.begin, self.end = [int(i) for i in interval.split('-')]
        self.size = self.end - self.begin
        self.string = "{}_{}_{}".format(self.chromosome, self.begin, self.end)


def make_windows(full_window, max_size=1e6):
    full_reg = Region(full_window)
    if full_reg.size > max_size:
        chunks = ceil(full_reg.size / max_size)
        chunksize = ceil(full_reg.size / chunks)
        return [
            Region("{}:{}-{}".format(
                full_reg.chromosome,
                full_reg.begin + i * chunksize,
                full_reg.begin + (i + 1) * chunksize))
            for i in range(chunks)]
    else:
        return [full_reg]


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
    parser.add_argument("--simplify",
                        help="simplify annotation track to show genes rather than transcripts",
                        action="store_true")
    parser.add_argument("--split",
                        help="split, rather than overlay the methylation tracks",
                        action="store_true")
    parser.add_argument("--smooth",
                        help="Smoothen the datapoints, but reduce the details",
                        type=int,
                        default=5)
    args = parser.parse_args()
    if not len(args.names) == len(args.methylation):
        sys.exit("INPUT ERROR: Expecting the same number of names as datasets!")
    return args
