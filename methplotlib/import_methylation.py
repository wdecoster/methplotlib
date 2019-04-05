import pandas as pd
import numpy as np


def read_meth_freq(filename, window):
    table = pd.read_csv(filename, sep="\t")
    table = table.loc[(table.chromosome == window.chromosome) &
                      (table.start > window.begin) &
                      (table.end < window.end)]
    table["pos"] = np.floor(table[['start', 'end']].mean(axis=1))
    return table.drop(columns=['start', 'end', 'num_motifs_in_group',
                               'called_sites', 'called_sites_methylated', 'group_sequence']) \
        .sort_values('pos') \
        .groupby('pos') \
        .mean() \
        .rolling(window=5, center=True) \
        .mean()


def get_data(methylation_files, window):
    """
    Import methylation frequency from all files in args.methylation
    within the window args.window
    """
    return [read_meth_freq(f, window) for f in methylation_files]
