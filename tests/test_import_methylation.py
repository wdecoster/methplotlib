import shutil

import pysam
import pytest

from methplotlib.import_methylation import parse_bedgraph, parse_bedmethyl, parse_nanopolish
from methplotlib.utils import Region

pytestmark = pytest.mark.skipif(shutil.which("tabix") is None, reason="tabix not on the PATH")

# `modkit pileup --only-tabs` records, as shared in issue #52
BEDMETHYL = """\
1\t10468\t10469\th\t1\t+\t10468\t10469\t255,0,0\t1\t0.00\t0\t0\t1\t1\t0\t0\t0
1\t10468\t10469\tm\t1\t+\t10468\t10469\t255,0,0\t1\t100.00\t1\t0\t0\t1\t0\t0\t0
1\t10469\t10470\th\t2\t-\t10469\t10470\t255,0,0\t2\t0.00\t0\t1\t1\t0\t0\t0\t0
1\t10469\t10470\tm\t2\t-\t10469\t10470\t255,0,0\t2\t50.00\t1\t1\t0\t0\t0\t0\t0
"""

BEDGRAPH = """\
1\t10468\t10469\t0.5
1\t10469\t10470\t1.0
"""

NANOPOLISH_FREQ = """\
chromosome\tstart\tend\tnum_motifs_in_group\tcalled_sites\tcalled_sites_methylated\tmethylated_frequency\tgroup_sequence
1\t10468\t10468\t1\t10\t5\t0.5\tACGTA
1\t10470\t10470\t1\t10\t10\t1.0\tACGTA
"""


def indexed(tmp_path, name, content, **index_args):
    path = tmp_path / name
    path.write_text(content)
    # bgzips the file and creates the .tbi index next to it
    return pysam.tabix_index(str(path), **index_args)


@pytest.fixture
def indexed_bedmethyl(tmp_path):
    return indexed(tmp_path, "test.bed", BEDMETHYL, preset="bed")


def corrupt_index(filename):
    with open(filename + ".tbi", "w") as index:
        index.write("not an index")


def test_indexed_bedmethyl_with_records(indexed_bedmethyl):
    mods = parse_bedmethyl(indexed_bedmethyl, "test", Region("1:10000-11000"))
    assert sorted(m.name for m in mods) == ["test_h", "test_m"]


@pytest.mark.parametrize("window", ["chrX:146,991,437-147,034,645", "1:500000-600000"])
def test_indexed_bedmethyl_without_records(indexed_bedmethyl, window):
    assert parse_bedmethyl(indexed_bedmethyl, "test", Region(window)) == []


@pytest.mark.parametrize(
    "name, content, index_args, parse",
    [
        ("test.bed", BEDMETHYL, {"preset": "bed"}, parse_bedmethyl),
        ("test.bedgraph", BEDGRAPH, {"preset": "bed"}, parse_bedgraph),
        (
            "freq.tsv",
            NANOPOLISH_FREQ,
            {"seq_col": 0, "start_col": 1, "end_col": 2, "line_skip": 1},
            lambda f, n, w: parse_nanopolish(f, "nanopolish_freq", n, w),
        ),
    ],
    ids=["bedmethyl", "bedgraph", "nanopolish_freq"],
)
def test_tabix_failure(tmp_path, name, content, index_args, parse):
    # a corrupt index makes tabix exit with an error and no output,
    # which should not be reported as a window without records
    filename = indexed(tmp_path, name, content, **index_args)
    corrupt_index(filename)
    with pytest.raises(SystemExit, match="tabix failed"):
        parse(filename, "test", Region("1:10000-11000"))
