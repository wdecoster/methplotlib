import pandas as pd
import pyranges as pr
import numpy as np
import sys
import logging
from methplotlib.utils import file_sniffer
import pysam


class Methylation(object):
    def __init__(self, table, data_type, name, called_sites):
        self.table = table
        self.data_type = data_type
        self.name = name
        self.called_sites = called_sites


def read_meth(filename, name, window, smoothen=5):
    """
    converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'
    smoothening the result by a rolling average

    input can also be raw data per read, optionally phased
    which will return a dataframe with 'read', 'chromosome', 'pos', 'log_lik_ratio', 'strand'
    """
    file_type = file_sniffer(filename)
    logging.info("File is of type {}".format(file_type))
    try:
        if file_type.startswith("nanopolish"):
            return parse_nanopolish(filename, file_type, name, window, smoothen=smoothen)
        elif file_type == "nanocompore":
            return parse_nanocompore(filename, name, window)
        elif file_type == "ont-cram":
            return parse_ont_cram(filename, name, window)
    except Exception:
        sys.stderr.write("\n\n\nInput file {} not recognized!\n".format(filename))
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise


def parse_nanopolish(filename, file_type, name, window, smoothen=5):
    table = pd.read_csv(filename, sep="\t")
    gr = pr.PyRanges(table.rename(columns={"start": "Start", "chromosome": "Chromosome",
                                           "end": "End", "Strand": "strand"}))
    logging.info("Read the file in a dataframe.")

    if window:
        gr = gr[window.chromosome, window.begin:window.end]
    gr.pos = np.floor(gr.drop().df[["Start", "End"]].mean(axis=1))
    table = gr.df

    if file_type in ['nanopolish_call', 'nanopolish_phased']:

        table = table.drop(columns=['Start', 'End', 'log_lik_methylated',
                                    'log_lik_unmethylated', 'num_calling_strands',
                                    'num_motifs', 'sequence'])
        return Methylation(
            table=table.sort_values(['read_name', 'pos']),
            data_type=file_type,
            name=name,
            called_sites=len(table))
    if file_type == "nanopolish_freq":
        called_sites = table.called_sites
        table = table.drop(columns=['Start', 'End', 'num_motifs_in_group',
                                    'called_sites', 'called_sites_methylated',
                                    'group_sequence'])
        return Methylation(
            table=table.sort_values('pos')
                       .groupby('pos')
                       .mean()
                       .rolling(window=smoothen, center=True)
                       .mean(),
            data_type=file_type,
            name=name,
            called_sites=called_sites.sum())


def parse_nanocompore(filename, name, window):
    def nanocompore_columns_of_interest(column):
        if column in ['pos', 'ref_id']:
            return True
        elif column.endswith('pvalue_context_2') or column.endswith('pvalue'):
            return True
        else:
            return False
    table = pd.read_csv(filename, sep="\t", usecols=nanocompore_columns_of_interest)
    if window:
        table = table[table["ref_id"] == window.chromosome]
    return Methylation(
        table=table.sort_values('pos')
                   .append({'pos': window.end}, ignore_index=True)
                   .drop(columns="ref_id")
                   .fillna(1.0),
        data_type='nanocompore',
        name=name,
        called_sites=len(table))


def parse_ont_cram(filename, name, window):
    cram = pysam.AlignmentFile(filename, "rc")
    data = []
    for read in cram.fetch(reference=window.chromosome, start=window.begin, end=window.end):
        if not read.is_supplementary and not read.is_secondary:
            mod, positions, quals = get_modified_reference_positions(read)
            for pos, qual in zip(positions, quals):
                if pos is not None:
                    data.append((read.query_name,
                                 '-' if read.is_reverse else '+',
                                 pos,
                                 qual,
                                 mod))
    return Methylation(
        table=pd.DataFrame(data, columns=['read_name', 'strand', 'pos', 'quality', 'mod'])
                .astype(dtype={'mod': 'category', 'quality': 'float'})
                .sort_values(['read_name', 'pos']),
        data_type="ont-cram",
        name=name,
        called_sites=len(data))


def get_modified_reference_positions(read):
    if read.has_tag('MM'):
        basemod = read.get_tag('MM').split(',', 1)[0]
        if '-' in basemod:
            sys.exit("ERROR: modifications on negative strand currently unsupported.")
        base, mod = basemod.split('+')
        deltas = [int(i) for i in read.get_tag('MM').split(',')[1:]]
        probabilities = phred_to_probability(read.get_tag('MP'))
        locations = np.cumsum(deltas) + np.concatenate(
            (np.zeros(shape=1),
             np.ones(shape=len(deltas) - 1))).astype('int')
        base_index = np.array(
            [i for i, letter in enumerate(read.get_forward_sequence()) if letter == base]
        )
        modified_bases = base_index[locations]
        refpos = np.array(read.get_reference_positions(full_length=True))
        if read.is_reverse:
            refpos = np.flipud(refpos)
            probabilities = probabilities[::-1]
        return (basemod, refpos[modified_bases], probabilities)
    else:
        return (None, [None], [None])


def errs_tab(n):
    """Generate list of error rates for qualities less than equal than n."""
    return [10**(q / -10) for q in range(n + 1)]


def phred_to_probability(quals, tab=errs_tab(128)):
    return [tab[ord(q) - 33] for q in quals]


def get_data(methylation_files, names, window, smoothen=5):
    """
    Import methylation data from all files in the list methylation_files

    Data can be either frequency or raw.

    data is extracted within the window args.window
    Frequencies are smoothened using a sliding window
    """
    return [read_meth(f, n, window, smoothen) for f, n in zip(methylation_files, names)]
