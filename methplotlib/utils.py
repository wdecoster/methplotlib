from argparse import ArgumentParser, SUPPRESS
import sys
from math import ceil
from methplotlib.version import __version__
from datetime import datetime as dt
from time import time
import logging
import plotly
from pathlib import Path
from itertools import chain
from plotly import subplots
from plotly.subplots import make_subplots


class Region(object):
    def __init__(self, region, fasta=None):
        if ":" in region:
            try:
                # the chromosome is deliberately kept as a string,
                # as that is how pyranges and tabix identify chromosomes
                self.chromosome, interval = region.replace(",", "").split(":")
                self.begin, self.end = [int(i) for i in interval.split("-")]
                self.start = self.begin
            except ValueError:
                sys.exit(
                    "\n\nERROR: Window (-w/--window) inproperly formatted, "
                    "examples of accepted formats are:\n"
                    "'chr5:150200605-150423790' or 'ENST00000647408'\n\n"
                )
            self.size = self.end - self.begin
            if not self.size > 0:
                sys.exit(
                    "\n\nERROR: Window (-w/--window) inproperly formatted, "
                    "begin of the interval has to be smaller than end\n\n"
                )
            self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        else:  # When region is an entire chromosome, contig or transcript
            if fasta is None:
                sys.exit(
                    "A fasta reference file is required if --window "
                    "is an entire chromosome, contig or transcript"
                )
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
            for i in range(chunks)
        ]
    else:
        return [reg]


def flatten(nested_list):
    return list(chain.from_iterable(nested_list))


def get_args():
    parser = ArgumentParser(description="plotting nanopolish methylation calls or frequency")
    parser.add_argument(
        "-v",
        "--version",
        help="Print version and exit.",
        action="version",
        version=f"methplotlib {__version__}",
    )
    parser.add_argument(
        "-m",
        "--methylation",
        nargs="+",
        help="data in nanopolish, nanocompore, ont-cram or bedgraph format",
        required=True if "--example" not in sys.argv else False,
    )
    parser.add_argument(
        "-n",
        "--names",
        nargs="+",
        help="names of datasets in --methylation",
        required=True if "--example" not in sys.argv else False,
    )
    parser.add_argument(
        "-w",
        "--window",
        help="window (region) to which the visualisation has to be restricted",
        required=True if "--example" not in sys.argv else False,
    )
    parser.add_argument("-g", "--gtf", help="add annotation based on a gtf file")
    parser.add_argument("-b", "--bed", help="add annotation based on a bed file")
    parser.add_argument(
        "-f",
        "--fasta",
        help="required when --window is an entire chromosome, contig or transcript",
    )
    parser.add_argument(
        "--simplify",
        help="simplify annotation track to show genes rather than transcripts",
        action="store_true",
    )
    parser.add_argument(
        "--split",
        help="split, rather than overlay the methylation tracks",
        action="store_true",
    )
    parser.add_argument("--static", help="Make a static image of the browser window")
    parser.add_argument(
        "--binary",
        help="Make the nanopolish plot ignorning log likelihood nuances",
        action="store_true",
    )
    parser.add_argument(
        "--smooth",
        help="Rolling window size for averaging frequency values",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--mods",
        help="Comma separated list of modifications to restrict to",
        default=None,
    )

    parser.add_argument(
        "--dotsize",
        help="Control the size of dots in the per read plots",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--minqual",
        help="The minimal phred quality to show [for bam/cram input only]",
        type=int,
        default=20,
    )
    parser.add_argument("--example", action="store_true", help="Show example command and exit.")
    parser.add_argument(
        "-o",
        "--outfile",
        help="File to write results to. "
        "Default: methylation_browser_{chr}_{start}_{end}.html. "
        "Use {region} as a shorthand for {chr}_{start}_{end} in the filename. "
        "Missing paths will be created.",
    )
    parser.add_argument(
        "-q",
        "--qcfile",
        help="File to write the qc report to. "
        "Default: The path in outfile prefixed with qc_, "
        "default is qc_report_methylation_browser_{chr}_{start}_{end}.html. "
        "Use {region} as a shorthand for {chr}_{start}_{end} in the filename. "
        "Missing paths will be created.",
    )
    parser.add_argument("--store", help=SUPPRESS, action="store_true")
    args = parser.parse_args()
    if not args.example and not len(args.names) == len(args.methylation):
        sys.exit("INPUT ERROR: Expecting the same number of names as datasets!")
    return args


def init_logs(args):
    """Initiate log file and log arguments."""
    start_time = dt.fromtimestamp(time()).strftime("%Y%m%d_%H%M")
    logname = "methplotlib_" + start_time + ".log"
    handlers = [logging.FileHandler(logname)]
    logging.basicConfig(format="%(asctime)s %(message)s", handlers=handlers, level=logging.INFO)
    logging.info(f"methplotlib {__version__} started.")
    py_version = sys.version.replace("\n", " ")
    logging.info(f"Python version is: {py_version}")
    logging.info(f"Arguments are: {args}")


def print_example():
    """
    Print an example command

    The example files are not shipped with the package, to keep it small,
    so the paths are relative to a clone of the repository
    """
    example = """
methplotlib -m examples/NA19240-methylation_ACTB_calls.tsv.gz \\
               examples/NA19240-methylation_ACTB_frequency.tsv.gz \\
            -n calls frequencies \\
            -w chr7:5,525,542-5,543,028 \\
            -g examples/GRCh38-ACTB-locus.gtf.gz \\
            --simplify \\
            -b examples/DNase_cluster_ACTB.bed.gz \\
            -o '{region}/example.html'""".strip()
    sys.stderr.write(
        "Example command, using the files in the examples folder of\n"
        "https://github.com/wdecoster/methplotlib\n\n"
    )
    print(example)
    sys.exit(0)


def is_gz_file(filepath):
    import binascii

    with open(filepath, "rb") as test_f:
        return binascii.hexlify(test_f.read(2)) == b"1f8b"


def is_cram_file(filepath):
    with open(filepath, "rb") as test_f:
        return test_f.read(4) == b"CRAM"


def is_bam_file(filepath):
    import gzip

    try:
        with gzip.open(filepath) as test_f:
            return test_f.read(3) == b"BAM"
    except OSError:
        return False


def open_text(filepath):
    """Open a file for reading text, transparently handling gzip compression."""
    if is_gz_file(filepath):
        import gzip

        return gzip.open(filepath, "rt")
    else:
        return open(filepath, "r")


def first_data_line(filename):
    """
    Return the first line of a file that contains data,
    skipping empty lines, comments and browser/track lines
    """
    with open_text(filename) as handle:
        for line in handle:
            if line.strip() and not line.startswith(("#", "track", "browser")):
                return line.rstrip("\n").rstrip("\r")
    sys.exit(f"\n\nERROR: File {filename} does not contain any data!\n")


def is_bedmethyl(fields):
    """
    Check if the fields of a line look like bedMethyl,
    which is BED9 with a modification code in column 4, followed by counts
    """
    return (
        len(fields) > 9
        and fields[1].isdigit()
        and fields[2].isdigit()
        and fields[5] in ["+", "-", "."]
    )


def is_bedgraph(fields):
    """
    Check if the fields of a line look like bedgraph:
    chromosome, start, end and a numeric value.
    A fifth column (coverage, as written by `modkit pileup --bedgraph`) is tolerated.
    """
    if len(fields) not in [4, 5] or not (fields[1].isdigit() and fields[2].isdigit()):
        return False
    try:
        float(fields[3])
        return True
    except ValueError:
        return False


def bedmethyl_separator(filename):
    """
    `modkit pileup` writes the first 9 (BED) columns tab separated
    and the remaining columns space separated, unless --only-tabs was used.
    Return the separator to read the file with.
    """
    line = first_data_line(filename)
    return "\t" if len(line.split("\t")) == len(line.split()) else r"\s+"


def unrecognized_file_error(filename, line):
    """Compose an error message explaining why a file could not be used."""
    fields = line.split("\t")
    whitespace_fields = line.split()
    return (
        f"\n\nERROR: Input file {filename} was not recognized!\n\n"
        f"Its first data line has {len(fields)} tab separated fields "
        f"({len(whitespace_fields)} when splitting on any whitespace):\n"
        f"{line[:200]}\n\n"
        "methplotlib accepts:\n"
        " - bedgraph: 4 columns (chromosome, start, end, value), "
        "or the 5 columns of `modkit pileup --bedgraph`\n"
        " - bedMethyl of `modkit pileup`: 18 columns, with or without --only-tabs\n"
        " - bedMethyl of `modbam2bed --extended`: 14 columns\n"
        " - bam or cram with MM/ML tags\n"
        " - nanopolish call-methylation output, optionally phased\n"
        " - the output of nanopolish/calculate_methylation_frequency.py\n"
        " - nanocompore results\n\n"
        "Frequent causes are a header line, a file that is comma or space separated "
        "throughout, or extra columns added by other tools.\n"
        "If your file should be supported, please open an issue at "
        "https://github.com/wdecoster/methplotlib/issues and share these lines.\n"
    )


def file_sniffer(filename):
    """
    Takes in a filename and tries to guess the input file type
    """
    if not Path(filename).is_file():
        sys.exit(f"\n\nERROR: File {filename} does not exist, please check the path!\n")
    if is_cram_file(filename):
        return "cram"
    if is_bam_file(filename):
        return "bam"
    header = first_data_line(filename)

    if "GMM_anova_pvalue" in header:
        return "nanocompore"
    if "log_lik_methylated" in header:
        if "PS" in header:
            return "nanopolish_phased"
        else:
            return "nanopolish_call"
    if "num_motifs_in_group" in header:
        return "nanopolish_freq"
    # bedMethyl files have no header and are tab separated,
    # except for modkit which separates the columns after BED9 by spaces
    # unless it was run with --only-tabs
    whitespace_fields = header.split()
    if is_bedmethyl(whitespace_fields):
        if len(whitespace_fields) == 14:  # file generated by `modbam2bed --extended`
            return "bedmethyl_extended"
        if len(whitespace_fields) == 18:  # file generated by `modkit pileup`
            return "bedmethyl"
    if is_bedgraph(header.split("\t")):
        return "bedgraph"
    sys.exit(unrecognized_file_error(filename, header))


def create_subplots(num_methrows, split, names=None, annotation=True):
    """
    Prepare the panels (rows * 1 column) for the subplots.
    If splitting: one row for each dataset, taking 90%/len(datasets) for heights
    If not: one row spanning 4 rows and taking 90% of the heights
    if annotation is True (bed or gtf) then add a row with height 10%
    """
    if split:
        return plotly.subplots.make_subplots(
            rows=num_methrows + annotation,
            cols=1,
            shared_xaxes=True,
            specs=[[{}] for i in range(num_methrows + annotation)],
            print_grid=False,
            subplot_titles=names,
            vertical_spacing=0.1 if num_methrows < 10 else 0.01,
            row_heights=[0.9 / num_methrows] * num_methrows + [0.1] * annotation,
        )
    else:
        return plotly.subplots.make_subplots(
            rows=num_methrows + annotation,
            cols=1,
            shared_xaxes=True,
            specs=[[{"rowspan": num_methrows}], [None], [None], [None]] + [[{}]] * annotation,
            print_grid=False,
            vertical_spacing=0.1 if num_methrows < 10 else 0.01,
            row_heights=[0.9, 0, 0, 0] + [0.1] * annotation,
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
        output.write(
            plotly.offline.plot(fig, output_type="div", show_link=False, include_plotlyjs="cdn")
        )
