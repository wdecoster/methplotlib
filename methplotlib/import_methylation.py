import pandas as pd
import numpy as np
import sys


def read_meth_freq(filename, window, smoothen=5):
    """
    converts a file from calculate_methylation_frequency to a pandas dataframe
    containing 'chromosome', 'pos', 'methylated_frequency'
    smoothening the result by a rolling average
    """
    table = pd.read_csv(filename, sep="\t")
    table = table.loc[(table.chromosome == window.chromosome) &
                      (table.start > window.begin) &
                      (table.end < window.end)]
    table["pos"] = np.floor(table[['start', 'end']].mean(axis=1))
    try:
        return table.drop(columns=['start', 'end', 'num_motifs_in_group',
                                   'called_sites', 'called_sites_methylated', 'group_sequence']) \
            .sort_values('pos') \
            .groupby('pos') \
            .mean() \
            .rolling(window=smoothen, center=True) \
            .mean()
    except Exception:
        sys.stderr.write("ERROR parsing {}\n\n\nDetailed error:\n".format(filename))
        raise


def get_data(methylation_files, window, smoothen):
    """
    Import methylation frequency from all files in the list methylation_files
    within the window args.window
    passing a smoothen parameter
    """
    return [read_meth_freq(f, window, smoothen) for f in methylation_files]
