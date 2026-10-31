import random

import numpy as np
import pysam

from methplotlib.import_methylation import reference_positions


def make_read(header, cigartuples, reference_start, reverse):
    read = pysam.AlignedSegment(header)
    read.query_name = "read"
    read.reference_id = 0
    read.reference_start = reference_start
    read.flag = 16 if reverse else 0
    query_length = sum(length for op, length in cigartuples if op in (0, 1, 4, 7, 8))
    read.query_sequence = "A" * query_length
    read.cigartuples = cigartuples
    return read


def test_reference_positions_match_pysam():
    """With every CIGAR operation, the same as get_reference_positions(full_length=True)"""
    rng = random.Random(1)
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 10**8}]})
    for _ in range(200):
        cigar = []
        if rng.random() < 0.5:
            cigar.append((5, rng.randrange(1, 50)))  # hard clip
        if rng.random() < 0.5:
            cigar.append((4, rng.randrange(1, 50)))  # soft clip
        for _ in range(rng.randrange(1, 100)):
            cigar.append((rng.choice([0, 7, 8]), rng.randrange(1, 40)))  # M, =, X
            cigar.append((rng.choice([1, 2, 3, 6]), rng.randrange(1, 6)))  # I, D, N, P
        cigar.append((0, 10))
        if rng.random() < 0.5:
            cigar.append((4, rng.randrange(1, 50)))
        if rng.random() < 0.5:
            cigar.append((5, rng.randrange(1, 50)))
        read = make_read(header, cigar, rng.randrange(0, 10**6), reverse=rng.random() < 0.5)
        expected = np.array(read.get_reference_positions(full_length=True), dtype=float)
        np.testing.assert_array_equal(reference_positions(read), expected)
