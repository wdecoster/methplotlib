import pandas as pd
import pyranges as pr
import numpy as np
import sys
import logging


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
    try:
        table = pd.read_csv(filename, sep="\t")
        gr = pr.PyRanges(table.rename(columns={"start": "Start", "chromosome": "Chromosome",
                                               "end": "End", "Strand": "strand"}))
        logging.info("Read the file in a dataframe.")

        if window:
            gr = gr[window.chromosome, window.begin:window.end]
        gr.pos = np.floor(gr.drop().df[["Start", "End"]].mean(axis=1))
        table = gr.df

        if 'log_lik_ratio' in table:  # indicating the file is 'raw' or 'phased'

            table = table.drop(columns=['Start', 'End', 'log_lik_methylated',
                                        'log_lik_unmethylated', 'num_calling_strands',
                                        'num_motifs', 'sequence'])
            if 'PS' in gr:  # indicating the file contains phased calls
                data_type = 'phased'
                logging.info("File contains phased raw data.")
            else:
                data_type = 'raw'
                logging.info("File contains raw data.")
            return Methylation(
                table=table.sort_values(['read_name', 'pos']),
                data_type=data_type,
                name=name,
                called_sites=len(table))
        else:  # assuming the file is from calculate_methylation_frequency
            logging.info("File contains frequency data.")
            called_sites = table.called_sites
            table = table.drop(columns=['Start', 'End', 'num_motifs_in_group',
                                        'called_sites', 'called_sites_methylated',
                                        'group_sequence'])
            return Methylation(
                table=table
                .sort_values('pos')
                .groupby('pos')
                .mean()
                .rolling(window=smoothen, center=True)
                .mean(),
                data_type="frequency",
                name=name,
                called_sites=called_sites.sum())
    except Exception:
        sys.stderr.write("\n\n\nInput file {} not recognized!\n".format(filename))
        sys.stderr.write("\n\n\nDetailed error:\n")
        raise


def get_data(methylation_files, names, window, smoothen=5):
    """
    Import methylation data from all files in the list methylation_files

    Data can be either frequency or raw.

    data is extracted within the window args.window
    Frequencies are smoothened using a sliding window
    """
    return [read_meth(f, n, window, smoothen) for f, n in zip(methylation_files, names)]
