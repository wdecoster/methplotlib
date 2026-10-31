import pandas as pd
import pyranges as pr
import itertools
import sys
from plotly.colors import DEFAULT_PLOTLY_COLORS as plcolors
import gzip
import logging
import os
from pathlib import Path

# the MANE Select transcripts of GRCh38, available as the built-in annotation "mane"
MANE_VERSION = "1.5"
MANE_FILE = f"MANE.GRCh38.v{MANE_VERSION}.ensembl_genomic.gtf.gz"
MANE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/"
    f"release_{MANE_VERSION}/{MANE_FILE}"
)


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


def cache_dir():
    """The directory to keep downloaded annotation in, following the XDG convention"""
    base = os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache"
    return Path(base) / "methplotlib"


def resolve_annotation(gtf):
    """
    Return the path of the annotation to use for --gtf

    "mane" is the built-in MANE Select annotation of GRCh38, downloaded once
    and kept in the cache directory, unless a file named "mane" exists.
    """
    if gtf is None or gtf.lower() != "mane" or Path(gtf).exists():
        return gtf
    path = cache_dir() / MANE_FILE
    if not path.is_file():
        import shutil
        import urllib.request

        sys.stderr.write(f"Downloading the MANE v{MANE_VERSION} annotation to {path}\n")
        path.parent.mkdir(parents=True, exist_ok=True)
        # downloading to a temporary name, so that an interrupted download is not used later
        partial = path.with_name(path.name + ".part")
        # the timeout is per connection attempt and read, not for the entire download,
        # and kept short as an unreachable IPv6 address is tried before IPv4
        try:
            with urllib.request.urlopen(MANE_URL, timeout=10) as response, open(partial, "wb") as out:
                shutil.copyfileobj(response, out)
        except OSError as e:
            partial.unlink(missing_ok=True)
            sys.exit(f"\n\nERROR: could not download the MANE annotation from {MANE_URL}:\n{e}\n")
        os.replace(partial, path)
    return str(path)


def open_gtf(gtff):
    """
    Open the gtf, using gzip if it's compressed
    based on extension
    """
    return gzip.open(gtff, "rt") if gtff.endswith(".gz") else open(gtff)


def good_record(fields, chromosome):
    """
    Filtering on the gtf lines
    by checking for right chromosome and right feature type

    Takes the first fields of a line, as splitting off the long attributes
    field for every line of an annotation is a waste.
    The chromosome has to match the entire field, as a prefix match would
    e.g. also accept chr10 up to chr19 and chr1_KI270706v1_random for chr1
    """
    return len(fields) > 2 and fields[0] == str(chromosome) and fields[2] in ["exon", "gene"]


def get_features(gtfline, type="gtf"):
    """
    Extract the desirable features from a gtf record
    """
    chromosome, _, feature, begin, end, _, strand, _, attributes = gtfline.split("\t")
    gene, transcript = parse_attributes(attributes.rstrip(), type=type)
    return [chromosome, feature, int(begin), int(end), strand, gene, transcript]


def parse_attributes(attributes, type="gtf"):
    """
    Parse the attributes string of gtf record
    Return the values corresponding to gene_name and transcript_id
    """
    attribute_delimiter = {"gtf": "; ", "gff": ";"}
    kv_delimiter = {"gtf": " ", "gff": "="}
    # the last attribute of a record keeps the trailing delimiter,
    # which would otherwise end up in the name shown when hovering
    info = {
        i.split(kv_delimiter[type])[0]: i.split(kv_delimiter[type])[1].replace('"', "").rstrip(";")
        for i in attributes.split(attribute_delimiter[type])
        if i.startswith(("gene_name", "transcript_id", "locus_tag"))
    }
    if "gene_name" in info.keys():
        return info.get("gene_name"), info.get("transcript_id")
    else:
        return info.get("locus_tag"), info.get("locus_tag")


def transcripts_in_window(df, window, feature="transcript"):
    """
    Return the transcript names of the exons overlapping the window

    Testing for overlap rather than for a begin or end within the window,
    as otherwise a feature spanning the entire window would be missed
    """
    return df.loc[
        (df["begin"] <= window.end) & (df["end"] >= window.begin),
        feature,
    ].unique()


def drop_genes_with_exons(df):
    """
    Remove the gene records of genes that also have exon records

    A gene record spans the entire gene, and would otherwise be drawn as one exon
    covering the introns. Genes without exon records keep their gene record.
    """
    genes_with_exons = df.loc[df["feature"] == "exon", "gene"].dropna().unique()
    return df[~((df["feature"] == "gene") & df["gene"].isin(genes_with_exons))]


def merge_intervals(intervals):
    """
    Merge overlapping and adjacent intervals of (begin, end) tuples,
    returned sorted by begin
    """
    merged = []
    for begin, end in sorted(intervals):
        if merged and begin <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((begin, end))
    return merged


def assign_colors_to_genes(transcripts):
    genes = dict.fromkeys([t.gene for t in transcripts])
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
    chromosomes_seen = set()
    records = []
    for line in open_gtf(gtff):
        if line.startswith("#") or not line.strip():
            continue
        # only the fields needed to decide, the attributes are split off
        # further down for the few lines that are of interest
        fields = line.split("\t", 3)
        if good_record(fields, window.chromosome):
            records.append(get_features(line, type=type))
        else:
            chromosomes_seen.add(fields[0])
    if not records and chromosomes_seen:
        sys.stderr.write(
            f"\nWARNING: no exons or genes of chromosome '{window.chromosome}' "
            f"in {gtff}.\n"
            f"The annotation contains {len(chromosomes_seen)} chromosomes, "
            f"such as {', '.join(sorted(chromosomes_seen)[:5])}.\n"
            "Do the chromosome names in --window match those in your annotation "
            "(e.g. '15' versus 'chr15')?\n"
        )
    df = pd.DataFrame(
        data=records,
        columns=["chromosome", "feature", "begin", "end", "strand", "gene", "transcript"],
    )
    logging.info("Loaded GTF file, processing...")
    df = drop_genes_with_exons(df)
    if simplify:
        res = []
        genes = transcripts_in_window(df, window, feature="gene")
        # grouping the features of interest once, rather than scanning the
        # entire frame again for every gene
        per_gene = dict(list(df[df["gene"].isin(genes)].groupby("gene")))
        for g in genes:
            gtable = per_gene.get(g)
            if gtable is not None and len(gtable):
                # a meta-transcript of the gene: the exons of all its transcripts, merged
                res.append(
                    Transcript(
                        transcript=gtable["gene"].tolist()[0],
                        gene=gtable["gene"].tolist()[0],
                        exon_tuples=merge_intervals(
                            gtable.loc[:, ["begin", "end"]].itertuples(index=False, name=None)
                        ),
                        strand=gtable["strand"].tolist()[0],
                    )
                )
        sys.stderr.write(f"Found {len(res)} gene(s) in the region.\n")
        logging.info(f"Found {len(res)} gene(s) in the region.\n")
    else:
        res = []
        transcripts = transcripts_in_window(df, window, feature="transcript")
        # grouping the features of interest once, rather than scanning the
        # entire frame again for every transcript
        per_transcript = dict(list(df[df["transcript"].isin(transcripts)].groupby("transcript")))
        for t in transcripts:
            tr = per_transcript.get(t)
            if tr is None or not len(tr):
                # an annotation without transcript_id, e.g. some gff flavors
                continue
            res.append(
                Transcript(
                    transcript=t,
                    gene=tr["gene"].tolist()[0],
                    exon_tuples=tr.loc[:, ["begin", "end"]]
                    .sort_values("begin")
                    .itertuples(index=False, name=None),
                    strand=tr["strand"].tolist()[0],
                )
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
    if annot_file.endswith((".gtf", ".gtf.gz")):
        return "gtf"
    elif annot_file.endswith((".gff", ".gff.gz")):
        return "gff"
    else:
        sys.exit(
            "ERROR: unrecognized extension of the annotation file.\n"
            "Supported are gtf, gtf.gz, gff and gff.gz"
        )


def parse_bed(bed, window):
    logging.info("Parsing BED file")
    gr = pr.read_bed(bed)[window.chromosome, window.begin : window.end]
    df = gr.unstrand().df
    df = df.drop(columns=["Chromosome", "Score", "Strand"], errors="ignore")
    if "Name" not in df.columns:
        df["Name"] = "noname"
    return df.itertuples(index=False, name=None)
