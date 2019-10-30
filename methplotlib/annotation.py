import pandas as pd
import itertools
import sys
from plotly.colors import DEFAULT_PLOTLY_COLORS as plcolors
import gzip
import logging


class Transcript(object):
    def __init__(self, transcript, gene, exon_tuples, strand):
        self.transcript = transcript
        self.gene = gene
        self.exon_tuples = list(exon_tuples)
        self.strand = strand
        self.marker = "triangle-right" if self.strand == "+" else "triangle-left"
        self.begin = min(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.end = max(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.color = ""


def open_gtf(gtff):
    """
    Open the gtf, using gzip if it's compressed
    based on extension
    """
    return gzip.open(gtff, 'rt') if gtff.endswith('.gz') else open(gtff)


def good_record(line, chromosome):
    '''
    Filtering on the gtf lines
    by checking for right chromosome and right feature type
    '''
    if line.startswith(chromosome) \
            and line.split('\t')[2] == 'exon':
        return True
    else:
        return False


def get_features(gtfline):
    """
    Extract the desirable features from a gtf record
    """
    chromosome, _, _, begin, end, _, strand, _, attributes = gtfline.split('\t')
    gene, transcript = parse_attributes(attributes.rstrip())
    return [chromosome, int(begin), int(end), strand, gene, transcript]


def parse_attributes(attributes):
    """
    Parse the attributes string of gtf record
    Return the values corresponding to gene_name and transcript_id
    """
    info = {i.split(' ')[0]: i.split(' ')[1].replace('"', '') for i in attributes.split(
        '; ') if i.startswith('gene_name') or i.startswith('transcript_id')}
    return info.get("gene_name"), info.get("transcript_id")


def transcripts_in_window(df, window, feature='transcript'):
    """
    Return the transcript names for which
    either the end or the begin of an exon is within the window
    """
    return df.loc[df['begin'].between(window.begin, window.end)
                  | df['end'].between(window.begin, window.end), feature] \
        .unique()


def parse_gtf(gtff, window, simplify=False):
    """
    Parse the gtff and select the relevant region as determined by the window
    return as Transcript objects
    """
    logging.info("Parsing GTF file...")
    df = pd.DataFrame(data=[get_features(line)
                            for line in open_gtf(gtff) if good_record(line, window.chromosome)],
                      columns=['chromosome', 'begin', 'end', 'strand', 'gene', 'transcript'])
    logging.info("Loaded GTF file, processing...")
    if simplify:
        df.drop_duplicates(subset=['chromosome', 'begin', 'end', 'gene'], inplace=True)
        res = []
        for g in transcripts_in_window(df, window, feature='gene'):
            gtable = df.loc[df["gene"] == g]
            if len(gtable):
                res.append(
                    Transcript(transcript=gtable["gene"].tolist()[0],
                               gene=gtable["gene"].tolist()[0],
                               exon_tuples=gtable.loc[:, ["begin", "end"]].sort_values("begin")
                               .itertuples(index=False,
                                           name=None),
                               strand=gtable["strand"].tolist()[0])
                )
        sys.stderr.write("Found {} gene(s) in the region.\n".format(len(res)))
        logging.info("Found {} gene(s) in the region.\n".format(len(res)))
    else:
        res = []
        for t in transcripts_in_window(df, window, feature="transcript"):
            tr = df.loc[df["transcript"] == t]
            res.append(
                Transcript(transcript=t,
                           gene=tr["gene"].tolist()[0],
                           exon_tuples=tr.loc[:, ["begin", "end"]].sort_values("begin")
                                                                  .itertuples(index=False,
                                                                              name=None),
                           strand=tr["strand"].tolist()[0])
            )
        sys.stderr.write("Found {} transcript(s) in the region.\n".format(len(res)))
        logging.info("Found {} transcript(s) in the region.\n".format(len(res)))
    assign_colors_to_genes(res)
    return res


def assign_colors_to_genes(transcripts):
    genes = set([t.gene for t in transcripts])
    colordict = {g: c for g, c in zip(genes, plcolors * 100)}
    for t in transcripts:
        t.color = colordict[t.gene]


def parse_bed(bed, window):
    logging.info("Parsing BED file")
    df = pd.read_csv(bed, sep="\t", names=['chromosome', 'begin', 'end', 'name', 'score', 'strand'])
    return df.loc[df['begin'].between(window.begin, window.end)
                  | df['end'].between(window.begin, window.end)] \
        .drop(columns=['chromosome', 'score', 'strand']) \
        .itertuples(index=False, name=None)
