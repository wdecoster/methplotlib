from argparse import ArgumentParser
import sys
from math import ceil
from methplotlib.version import __version__
from datetime import datetime as dt
from time import time
import logging
import binascii
import gzip


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
                        required=True if "--example" not in sys.argv else False)
    parser.add_argument("-n", "--names",
                        nargs='+',
                        help="names of datasets in --methylation",
                        required=True if "--example" not in sys.argv else False)
    parser.add_argument("-w", "--window",
                        help="window (region) to which the visualisation has to be restricted",
                        required=True if "--example" not in sys.argv else False)
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
    parser.add_argument("--dotsize",
                        help="Control the size of dots in the per read plots",
                        type=int,
                        default=4)
    parser.add_argument("--example",
                        action="store_true",
                        help="Show example command and exit.")
    parser.add_argument("-o", "--outfile",
                        help="File to write results to. "
                             "Default: methylation_browser_{chr}_{start}_{end}.html. "
                             "Use {region} as a shorthand for {chr}_{start}_{end} in the filename. "
                             "Missing paths will be created.")
    parser.add_argument("-q", "--qcfile",
                        help="File to write the qc report to. "
                             "Default: The path in outfile prefixed with qc_, "
                             "default is qc_report_methylation_browser_{chr}_{start}_{end}.html. "
                             "Use {region} as a shorthand for {chr}_{start}_{end} in the filename. "
                             "Missing paths will be created.")

    args = parser.parse_args()
    if not args.example and not len(args.names) == len(args.methylation):
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


def print_example():
    import pkg_resources
    meth = pkg_resources.resource_filename("methplotlib", "examples/ACTB_calls.tsv.gz")
    meth_freq = pkg_resources.resource_filename("methplotlib", "examples/meth_freq.tsv.gz")
    bed = pkg_resources.resource_filename("methplotlib", "examples/DNAse_cluster.bed.gz")
    annotation = pkg_resources.resource_filename("methplotlib", "examples/g38_locus.gtf.gz")

    example = """
methplotlib -m {meth} \\
               {meth_freq} \\
            -n calls frequencies \\
            -w chr7:5,525,542-5,543,028 \\
            -g {annotation} \\
            --simplify \\
            -b {bed} \\
            -o '{{region}}/example.html'""".strip().format(meth=meth, meth_freq=meth_freq,
                                                           annotation=annotation, bed=bed)
    print(example)
    sys.exit(0)


def is_gz_file(filepath):
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def file_sniffer(filename):
    """
    Takes in a filename and tries to guess the input file type
    """
    if is_gz_file(filename):
        header = gzip.open(filename, 'rt').readline()
    else:
        header = open(filename, 'r').readline()
    if "GMM_anova_pvalue" in header:
        return "nanocompore"
    if "log_lik_methylated" in header:
        if "PS" in header:
            return "nanopolish_phased"
        else:
            return "nanopolish_call"
    if "num_motifs_in_group" in header:
        return "nanopolish_freq"
    else:
        sys.exit("\n\n\nInput file {} not recognized!\n".format(filename))
