import pandas as pd
import pyranges as pr
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
            and line.split('\t')[2] in ['exon', 'gene']:
        return True
    else:
        return False


def get_features(gtfline, type="gtf"):
    """
    Extract the desirable features from a gtf record
    """
    chromosome, _, _, begin, end, _, strand, _, attributes = gtfline.split('\t')
    gene, transcript = parse_attributes(attributes.rstrip(), type=type)
    return [chromosome, int(begin), int(end), strand, gene, transcript]


def parse_attributes(attributes, type="gtf"):
    """
    Parse the attributes string of gtf record
    Return the values corresponding to gene_name and transcript_id
    """
    attribute_delimiter = {'gtf': '; ', 'gff': ';'}
    kv_delimiter = {'gtf': ' ', 'gff': '='}
    info = {i.split(kv_delimiter[type])[0]: i.split(kv_delimiter[type])[1].replace('"', '')
            for i in attributes.split(attribute_delimiter[type])
            if i.startswith(('gene_name', 'transcript_id', 'locus_tag'))}
    if "gene_name" in info.keys():
        return info.get("gene_name"), info.get("transcript_id")
    else:
        return info.get("locus_tag"), info.get("locus_tag")


def transcripts_in_window(df, window, feature='transcript'):
    """
    Return the transcript names for which
    either the end or the begin of an exon is within the window
    """
    return df.loc[df['begin'].between(window.begin, window.end)
                  | df['end'].between(window.begin, window.end), feature] \
        .unique()


def assign_colors_to_genes(transcripts):
    genes = set([t.gene for t in transcripts])
    colordict = {g: c for g, c in zip(genes, plcolors * 100)}
    for t in transcripts:
        t.color = colordict[t.gene]


def parse_annotation(gtff, window, simplify=False):
    """
    Parse the gtff and select the relevant region as determined by the window
    return as Transcript objects
    """
    type = annot_file_sniffer(gtff)
    logging.info(f"Parsing {type} file...")
    df = pd.DataFrame(data=[get_features(line, type=type)
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
        sys.stderr.write(f"Found {len(res)} gene(s) in the region.\n")
        logging.info(f"Found {len(res)} gene(s) in the region.\n")
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
        sys.stderr.write(f"Found {len(res)} transcript(s) in the region.\n")
        logging.info(f"Found {len(res)} transcript(s) in the region.\n")
    assign_colors_to_genes(res)
    return res


def annot_file_sniffer(annot_file):
    """
    Figure out type of annotation file

    Right not just lazily focus on the extension
    """
    if annot_file.endswith(('.gtf', '.gtf.gz')):
        return 'gtf'
    elif annot_file.endswith(('.gff', '.gff.gz')):
        return 'gff'
    else:
        sys.exit("ERROR: unrecognized extension of the annotation file.\n"
                 "Supported are gtf, gtf.gz, gff and gff.gz")


def parse_bed(bed, window):
    logging.info("Parsing BED file")
    gr = pr.read_bed(bed)[window.chromosome, window.begin:window.end]
    df = gr.unstrand().df
    df = df.drop(columns=["Chromosome", "Score", "Strand"], errors='ignore')
    if "Name" not in df.columns:
        df["Name"] = "noname"
    return df.itertuples(index=False, name=None)
