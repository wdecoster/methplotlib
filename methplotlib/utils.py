from argparse import ArgumentParser, SUPPRESS
import sys
from math import ceil
from methplotlib.version import __version__
from datetime import datetime as dt
from time import time
import logging
import plotly
from pathlib import Path


class Region(object):
    def __init__(self, region, fasta=None):
        if ':' in region:
            try:
                self.chromosome, interval = region.replace(',', '').split(':')
                self.begin, self.end = [int(i) for i in interval.split('-')]
            except ValueError:
                sys.exit("\n\nERROR: Window (-w/--window) inproperly formatted, "
                         "examples of accepted formats are:\n"
                         "'chr5:150200605-150423790' or 'ENST00000647408'\n\n")
            self.size = self.end - self.begin
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        else:  # When region is an entire chromosome, contig or transcript
            if fasta is None:
                sys.exit("A fasta reference file is required if --window "
                         "is an entire chromosome, contig or transcript")
            else:
                from pyfaidx import Fasta
                self.chromosome = region
                self.begin = 0
                self.string = region
                self.end = len(Fasta(fasta)[region])
                self.size = self.end


def make_windows(full_window, max_size=1e6, fasta=None):
    reg = Region(full_window, fasta)
    if reg.size > max_size:
        chunks = ceil(reg.size / max_size)
        chsize = ceil(reg.size / chunks)
        return [
            Region(f"{reg.chromosome}:{reg.begin + i * chsize}-{reg.begin + (i + 1) * chsize}")
            for i in range(chunks)]
    else:
        return [reg]


def get_args():
    parser = ArgumentParser(description="plotting nanopolish methylation calls or frequency")
    parser.add_argument("-v", "--version",
                        help="Print version and exit.",
                        action="version",
                        version=f'methplotlib {__version__}')
    parser.add_argument("-m", "--methylation",
                        nargs='+',
                        help="data in nanopolish, nanocompore, ont-cram or bedgraph format",
                        required=True if "--example" not in sys.argv else False)
    parser.add_argument("-n", "--names",
                        nargs='+',
                        help="names of datasets in --methylation",
                        required=True if "--example" not in sys.argv else False)
    parser.add_argument("-w", "--window",
                        help="window (region) to which the visualisation has to be restricted",
                        required=True if "--example" not in sys.argv else False)
    parser.add_argument("-g", "--gtf",
                        help="add annotation based on a gtf file")
    parser.add_argument("-b", "--bed",
                        help="add annotation based on a bed file")
    parser.add_argument("-f", "--fasta",
                        help="required when --window is an entire chromosome, contig or transcript")
    parser.add_argument("--simplify",
                        help="simplify annotation track to show genes rather than transcripts",
                        action="store_true")
    parser.add_argument("--split",
                        help="split, rather than overlay the methylation tracks",
                        action="store_true")
    parser.add_argument("--static",
                        help="Make a static image of the browser window")
    parser.add_argument("--binary",
                        help="Make the nanopolish plot ignorning log likelihood nuances",
                        action="store_true")
    parser.add_argument("--smooth",
                        help="Rolling window size for averaging frequency values",
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
    parser.add_argument("--store",
                        help=SUPPRESS,
                        action="store_true")
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
    logging.info(f'methplotlib {__version__} started.')
    py_version = sys.version.replace('\n', ' ')
    logging.info(f'Python version is: {py_version}')
    logging.info(f'Arguments are: {args}')


def print_example():
    import pkg_resources
    meth = pkg_resources.resource_filename("methplotlib", "examples/ACTB_calls.tsv.gz")
    meth_freq = pkg_resources.resource_filename("methplotlib", "examples/meth_freq.tsv.gz")
    bed = pkg_resources.resource_filename("methplotlib", "examples/DNAse_cluster.bed.gz")
    annotation = pkg_resources.resource_filename("methplotlib", "examples/g38_locus.gtf.gz")

    example = f"""
methplotlib -m {meth} \\
               {meth_freq} \\
            -n calls frequencies \\
            -w chr7:5,525,542-5,543,028 \\
            -g {annotation} \\
            --simplify \\
            -b {bed} \\
            -o '{{region}}/example.html'""".strip()
    print(example)
    sys.exit(0)


def is_gz_file(filepath):
    import binascii
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'


def is_cram_file(filepath):
    with open(filepath, 'rb') as test_f:
        return test_f.read(4) == b'CRAM'


def file_sniffer(filename):
    """
    Takes in a filename and tries to guess the input file type
    """
    if not Path(filename).is_file():
        sys.exit(f"\n\nERROR: File {filename} does not exist, please check the path!\n")
    if is_cram_file(filename):
        return "ont-cram"
    if is_gz_file(filename):
        import gzip
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
    if len(header.split('\t')) == 4:
        return 'bedgraph'  # there is no header, but the file is tab separated and has 4 fields
    sys.exit(f"\n\n\nInput file {filename} not recognized!\n")


def create_subplots(num_methrows, split, names=None, annotation=True):
    '''
    Prepare the panels (rows * 1 column) for the subplots.
    If splitting: one row for each dataset, taking 90%/len(datasets) for heights
    If not: one row spanning 4 rows and taking 90% of the heights
    if annotation is True (bed or gtf) then add a row with height 10%
    '''
    if split:
        return plotly.subplots.make_subplots(
            rows=num_methrows + annotation,
            cols=1,
            shared_xaxes=True,
            specs=[[{}] for i in range(num_methrows + annotation)],
            print_grid=False,
            subplot_titles=names,
            vertical_spacing=0.1 if num_methrows < 10 else 0.01,
            row_heights=[0.9 / num_methrows] * num_methrows + [0.1] * annotation

        )
    else:
        return plotly.subplots.make_subplots(
            rows=num_methrows + annotation,
            cols=1,
            shared_xaxes=True,
            specs=[[{'rowspan': num_methrows}], [None], [None], [None]] + [[{}]] * annotation,
            print_grid=False,
            vertical_spacing=0.1 if num_methrows < 10 else 0.01,
            row_heights=[0.9, 0, 0, 0] + [0.1] * annotation
        )


def create_browser_output(fig, outfile, window):
    if outfile is None:
        outfile = f"methylation_browser_{window.string}.html"
    else:
        from pathlib import Path

        outfile = outfile.format(region=window.string)
        p = Path(outfile)
        Path.mkdir(p.parent, exist_ok=True, parents=True)

    if outfile.endswith(".html"):
        write_html_output(fig, outfile)
    else:
        try:
            fig.write_image(outfile)
        except ValueError as e:
            sys.stderr.write("\n\nERROR: creating the image in this file format failed.\n")
            sys.stderr.write("ERROR: creating in default html format instead.\n")
            sys.stderr.write("ERROR: additional packages required. Detailed error:\n")
            sys.stderr.write(str(e))
            write_html_output(fig, outfile)


def write_html_output(fig, outfile):
    with open(outfile, "w+") as output:
        output.write(plotly.offline.plot(fig,
                                         output_type="div",
                                         show_link=False,
                                         include_plotlyjs='cdn'))
