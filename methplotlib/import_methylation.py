import pandas as pd
import pyranges as pr
import numpy as np
import sys
import logging
from methplotlib.utils import file_sniffer


class Methylation(object):
    def __init__(self, table, data_type, name, called_sites):
        self.table = table
        self.data_type = data_type
        self.name = name
        self.called_sites = called_sites


def get_data(methylation_files, names, window, smoothen=5):
    """
    Import methylation data from all files in the list methylation_files

    Data can in various formats
    - nanopolish raw, phased, frequency
    - cram
    - nanocompore
    - bedgraph

    data is extracted within the window args.window
    Frequencies are smoothened using a sliding window
    """
    return [read_meth(f, n, window, smoothen) for f, n in zip(methylation_files, names)]


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
    logging.info(f"File {filename} is of type {file_type}")
    try:
        if file_type.startswith("nanopolish"):
            return parse_nanopolish(filename, file_type, name, window, smoothen=smoothen)
        elif file_type == "nanocompore":
            return parse_nanocompore(filename, name, window)
        elif file_type == "ont-cram":
            return parse_ont_cram(filename, name, window)
        elif file_type == 'bedgraph':
            return parse_bedgraph(filename, name, window)
    except Exception as e:
        logging.error(f"Error processing {filename}.")
        logging.error(e, exc_info=True)
        sys.stderr.write(f"\n\n\nError processing {filename}!\n")
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise


def parse_nanopolish(filename, file_type, name, window, smoothen=5):
    if window:
        from pathlib import Path
        if Path(filename + '.tbi').is_file():
            import subprocess
            import gzip
            logging.info(f"Reading {filename} using a tabix stream.")
            region = f"{window.chromosome}:{window.begin}-{window.end}"
            try:
                tabix_stream = subprocess.Popen(['tabix', filename, region],
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE)
            except FileNotFoundError as e:
                logging.error("Error when opening a tabix stream.")
                logging.error(e, exc_info=True)
                sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
                sys.stderr.write("Is tabix installed and on the PATH?.")
                raise
            header = gzip.open(filename, 'rt').readline().rstrip().split('\t')
            table = pd.read_csv(tabix_stream.stdout, sep='\t', header=None, names=header)
        else:
            logging.info(f"Reading {filename} slowly by splitting the file in chunks.")
            sys.stderr.write(f"\nReading {filename} would be faster with bgzip and tabix.\n")
            if file_type in ['nanopolish_call', 'nanopolish_phased']:
                sys.stderr.write("Please index with 'tabix -S1 -s1 -b3 -e4'.\n")
            else:
                sys.stderr.write("Please index with 'tabix -S1 -s1 -b2 -e3'.\n")
            iter_csv = pd.read_csv(filename, sep="\t", iterator=True, chunksize=1e6)
            table = pd.concat([chunk[chunk['chromosome'] == window.chromosome]
                               for chunk in iter_csv])
    else:
        table = pd.read_csv(filename, sep="\t")
    gr = pr.PyRanges(table.rename(columns={"start": "Start", "chromosome": "Chromosome",
                                           "end": "End", "Strand": "strand"}))
    logging.info("Read the file in a dataframe.")

    if window:
        gr = gr[window.chromosome, window.begin:window.end]
    try:
        gr.pos = np.floor(gr.drop().df[["Start", "End"]].mean(axis=1))
    except KeyError:
        sys.stderr.write(f"\n\n\nProblem parsing nanopolish file {filename}!\n")
        sys.stderr.write("Could it be that there are no calls in your selected window?\n")
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise

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


def parse_bedgraph(filename, name, window):
    if window:
        from pathlib import Path
        if Path(filename + '.tbi').is_file():
            import subprocess
            logging.info(f"Reading {filename} using a tabix stream.")
            region = f"{window.chromosome}:{window.begin}-{window.end}"
            try:
                tabix_stream = subprocess.Popen(['tabix', filename, region],
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.PIPE)
            except FileNotFoundError as e:
                logging.error("Error when opening a tabix stream.")
                logging.error(e, exc_info=True)
                sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
                sys.stderr.write("Is tabix installed and on the PATH?.")
                raise
            table = pd.read_csv(tabix_stream.stdout, sep='\t', header=None,
                                names=['Chromosome', 'Start', 'End', 'Value'])
        else:
            logging.info(f"Reading {filename} slowly by splitting the file in chunks.")
            sys.stderr.write(
                f"\nReading {filename} would be faster with bgzip and 'tabix -p bed'.\n")
            iter_csv = pd.read_csv(filename, sep="\t", iterator=True, chunksize=1e6, header=None,
                                   names=['Chromosome', 'Start', 'End', 'Value'])
            table = pd.concat([chunk[chunk['Chromosome'] == window.chromosome]
                               for chunk in iter_csv])
    else:
        table = pd.read_csv(filename, sep="\t", header=None,
                            names=['Chromosome', 'Start', 'End', 'Value'])
    gr = pr.PyRanges(table)
    logging.info("Read the file in a dataframe.")
    if window:
        gr = gr[window.chromosome, window.begin:window.end]
        if len(gr.df) == 0:
            sys.exit(f"No records for {filename} in {window.string}!\n")
    return Methylation(
        table=gr.df.sort_values('Start'),
        data_type="bedgraph",
        name=name,
        called_sites=len(table))


def parse_ont_cram(filename, name, window):
    import pysam
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
