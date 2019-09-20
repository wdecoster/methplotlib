from argparse import ArgumentParser
import sys
from math import ceil
from methplotlib.version import __version__
from datetime import datetime as dt
from time import time
import logging


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
    parser = ArgumentParser(description="plotting nanopolish methylation calls or frequency")
    parser.add_argument("-v", "--version",
                        help="Print version and exit.",
                        action="version",
                        version='methplotlib {}'.format(__version__))
    parser.add_argument("-m", "--methylation",
                        nargs='+',
                        help="nanopolish methylation calls or frequency output",
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
    parser.add_argument("-b", "--bed",
                        help="add annotation based on a bed file matching to your reference genome")
    parser.add_argument("--simplify",
                        help="simplify annotation track to show genes rather than transcripts",
                        action="store_true")
    parser.add_argument("--split",
                        help="split, rather than overlay the methylation tracks",
                        action="store_true")
    parser.add_argument("--smooth",
                        help="When plotting frequencies points are averaged using a rolling window",
                        type=int,
                        default=5)
    args = parser.parse_args()
    if not len(args.names) == len(args.methylation):
        sys.exit("INPUT ERROR: Expecting the same number of names as datasets!")
    return args


def init_logs(args):
    """Initiate log file and log arguments."""
    start_time = dt.fromtimestamp(time()).strftime('%Y%m%d_%H%M')
    logname = "methplotlib_" + start_time + ".log"
    handlers = [logging.FileHandler(logname)]
    logging.basicConfig(
        format='%(asctime)s %(message)s',
        handlers=handlers,
        level=logging.INFO)
    logging.info('methplotlib {} started.\nPython version is: {}\nArguments are: {}'.format(
        __version__, sys.version.replace('\n', ' '), args))
