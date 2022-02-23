#! /usr/bin/env python
# Part of nanopolish
# Copied here for convenience
# Slightly edited by wdecoster


import sys
import csv
import argparse
import gzip


class SiteStats:
    def __init__(self, g_size, g_seq):
        self.num_reads = 0
        self.called_sites = 0
        self.called_sites_methylated = 0
        self.group_size = g_size
        self.sequence = g_seq


def update_call_stats(key, num_called_cpg_sites, is_methylated, sequence):
    if key not in sites:
        sites[key] = SiteStats(num_called_cpg_sites, sequence)

    sites[key].num_reads += 1
    sites[key].called_sites += num_called_cpg_sites
    if is_methylated > 0:
        sites[key].called_sites_methylated += num_called_cpg_sites


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate methylation frequency per site"
    )
    parser.add_argument("-c", "--call-threshold", type=float, required=False, default=2)
    parser.add_argument("-i", "--input", type=str, required=False)
    parser.add_argument("-s", "--split-groups", action="store_true")
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="file doesn't have the expected header line",
    )
    args = parser.parse_args()
    assert args.call_threshold is not None

    sites = dict()

    if args.input:
        if args.input.endswith(".gz"):
            in_fh = gzip.open(args.input, "rt")
        else:
            in_fh = open(args.input)
    else:
        in_fh = sys.stdin
    if args.no_header:
        reader = csv.DictReader(
            in_fh,
            delimiter="\t",
            fieldnames=[
                "chromosome",
                "strand",
                "start",
                "end",
                "read_name",
                "log_lik_ratio",
                "log_lik_methylated",
                "log_lik_unmethylated",
                "num_calling_strands",
                "num_motifs",
                "sequence",
                "PS",
                "HP",
            ],
        )
    else:
        reader = csv.DictReader(in_fh, delimiter="\t")
    for record in reader:

        num_sites = int(record["num_motifs"])
        llr = float(record["log_lik_ratio"])

        # Skip ambiguous call
        if abs(llr) < args.call_threshold * num_sites:
            continue
        sequence = record["sequence"]

        is_methylated = llr > 0

        # if this is a multi-cpg group and split_groups is set, break up these sites
        if args.split_groups and num_sites > 1:
            c = str(record["chromosome"])
            s = int(record["start"])
            e = int(record["end"])

            # find the position of the first CG dinucleotide
            sequence = record["sequence"]
            cg_pos = sequence.find("CG")
            first_cg_pos = cg_pos
            while cg_pos != -1:
                key = (c, s + cg_pos - first_cg_pos, s + cg_pos - first_cg_pos)
                update_call_stats(key, 1, is_methylated, "split-group")
                cg_pos = sequence.find("CG", cg_pos + 1)
        else:
            key = (str(record["chromosome"]), int(record["start"]), int(record["end"]))
            update_call_stats(key, num_sites, is_methylated, sequence)

    # header
    print(
        "\t".join(
            [
                "chromosome",
                "start",
                "end",
                "num_motifs_in_group",
                "called_sites",
                "called_sites_methylated",
                "methylated_frequency",
                "group_sequence",
            ]
        )
    )

    sorted_keys = sorted(sites.keys(), key=lambda x: x)
    if len(sorted_keys) == 0:
        sys.exit("ERROR: No sites found for calculating frequencies!")

    for key in sorted_keys:
        if sites[key].called_sites > 0:
            (c, s, e) = key
            f = float(sites[key].called_sites_methylated) / sites[key].called_sites
            print(
                "%s\t%s\t%s\t%d\t%d\t%d\t%.3f\t%s"
                % (
                    c,
                    s,
                    e,
                    sites[key].group_size,
                    sites[key].called_sites,
                    sites[key].called_sites_methylated,
                    f,
                    sites[key].sequence,
                )
            )
