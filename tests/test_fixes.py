import array
import gzip

import numpy as np
import pandas as pd
import pysam

from methplotlib.annotation import merge_intervals, parse_annotation
from methplotlib.import_methylation import get_modified_reference_positions, parse_bedgraph
from methplotlib.plots import make_per_read_meth_traces_phred
from methplotlib.utils import Region, make_windows


def make_read(sequence, mm, ml, reverse=False):
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 10**6}]})
    read = pysam.AlignedSegment(header)
    read.query_name = "read"
    read.reference_id = 0
    read.reference_start = 100
    read.flag = 16 if reverse else 0
    read.query_sequence = sequence
    read.cigartuples = [(0, len(sequence))]
    read.set_tag("MM", mm)
    read.set_tag("ML", array.array("B", ml))
    return read


def test_combined_codes_are_split():
    """C+mh has the likelihoods of m and h interleaved per position"""
    read = make_read("ACGTCGTCGA", "C+mh,0,1;A+a,0;", [10, 20, 30, 40, 50])
    blocks = {block[4]: block for block in get_modified_reference_positions(read)}
    assert sorted(blocks) == ["A+a", "C+h", "C+m"]
    assert list(blocks["C+m"][2]) == [101, 107]
    assert list(blocks["C+m"][3]) == [10, 30]
    assert list(blocks["C+h"][2]) == [101, 107]
    assert list(blocks["C+h"][3]) == [20, 40]
    assert list(blocks["A+a"][3]) == [50]


def test_combined_codes_keep_suffix_and_chebi():
    read = make_read("ACGTCGTCGA", "C+mh.,0;C+21839?,1;", [10, 20, 30])
    mods = [(block[4], list(block[3])) for block in get_modified_reference_positions(read)]
    assert mods == [("C+m.", [10]), ("C+h.", [20]), ("C+21839", [30])]


def test_reads_without_calls_keep_their_strand():
    calls = pd.DataFrame(
        {"read_name": ["fwd"], "strand": ["+"], "pos": [150], "quality": [200.0], "mod": ["C+m"]}
    )
    reads = pd.DataFrame(
        [("fwd", 100, 200, "+"), ("rev_without_calls", 120, 220, "-"), ("fwd_without_calls", 300, 400, "+")],
        columns=["read_name", "posmin", "posmax", "strand"],
    ).set_index("read_name")
    line = make_per_read_meth_traces_phred(calls, reads)[0]
    symbols = dict(zip(line.x[::3], line.marker.symbol[::3]))
    assert symbols == {100: "triangle-right", 120: "triangle-left", 300: "triangle-right"}


def write_gtf(path, records):
    with gzip.open(path, "wt") as out:
        for feature, begin, end, attributes in records:
            out.write(f"chr1\ttest\t{feature}\t{begin}\t{end}\t.\t+\t.\t{attributes}\n")


def test_simplify_merges_exons_and_ignores_gene_records(tmp_path):
    gtf = str(tmp_path / "test.gtf.gz")
    gene = 'gene_id "G1"; gene_name "GENE";'
    write_gtf(
        gtf,
        [
            ("gene", 100, 1000, gene),
            ("exon", 100, 200, gene + ' transcript_id "T1";'),
            ("exon", 500, 600, gene + ' transcript_id "T1";'),
            ("exon", 150, 250, gene + ' transcript_id "T2";'),
            ("exon", 900, 1000, gene + ' transcript_id "T2";'),
            ("gene", 2000, 3000, 'gene_id "G2"; gene_name "NOEXONS";'),
        ],
    )
    genes = parse_annotation(gtf, Region("chr1:1-5000"), simplify=True)
    assert [(g.gene, g.exon_tuples) for g in genes] == [
        ("GENE", [(100, 250), (500, 600), (900, 1000)]),
        ("NOEXONS", [(2000, 3000)]),
    ]
    transcripts = parse_annotation(gtf, Region("chr1:1-5000"), simplify=False)
    assert [(t.transcript, t.exon_tuples) for t in transcripts] == [
        ("T1", [(100, 200), (500, 600)]),
        ("T2", [(150, 250), (900, 1000)]),
    ]


def test_merge_intervals():
    assert merge_intervals([(5, 8), (1, 3), (4, 4), (10, 12), (11, 20)]) == [(1, 8), (10, 20)]


def test_header_lines_are_skipped(tmp_path):
    bedgraph = tmp_path / "test.bedgraph"
    bedgraph.write_text("track type=bedGraph\nchr1\t10\t20\t0.5\nchr1\t20\t30\t0.7\n")
    [result] = parse_bedgraph(str(bedgraph), "test", Region("chr1:1-100"))
    assert list(result.table["Value"]) == [0.5, 0.7]


def test_last_window_ends_at_the_end():
    windows = make_windows("chr1:0-2500001")
    assert windows[-1].end == 2500001
    assert all(a.end == b.begin for a, b in zip(windows, windows[1:]))
