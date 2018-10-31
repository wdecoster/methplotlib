import pandas as pd
import numpy as np
import sys
import plotly
import plotly.graph_objs as go
from argparse import ArgumentParser


def parse_region(region):
    chromosome, interval = region.split(':')
    begin, end = [int(i) for i in interval.split('-')]
    return chromosome, begin, end


def read_meth_freq(filename, window):
    table = pd.read_csv(filename, sep="\t")
    chromosome, begin, end = parse_region(window)
    table = table.loc[(table.chromosome == chromosome) & (table.start > begin) & (table.end < end)]
    table["pos"] = np.floor(table[['start', 'end']].mean(axis=1))
    table = table.drop(columns=['start', 'end', 'num_cpgs_in_group',
                                'called_sites', 'called_sites_methylated', 'group_sequence']) \
        .sort_values('pos') \
        .groupby('pos') \
        .mean() \
        .rolling(window=5, center=True) \
        .mean()
    return table


def get_data(args):
    return [read_meth_freq(f, args.window) for f in args.methylation]


def meth_browser(methlist, names):
    data = [go.Scatter(x=a.index, y=a["methylated_frequency"],
                       mode='lines+markers',
                       name=n) for a, n in zip(methlist, names)]

    html = plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title="Methylation")
         },
        output_type="div",
        show_link=False)
    with open("methylation_browser.html", 'w') as output:
        output.write(html)


def main():
    args = get_args()
    methlist = get_data(args)
    meth_browser(methlist, args.names)


def get_args():
    parser = ArgumentParser(description="plotting methylation for two alleles separately")
    parser.add_argument("-m", "--methylation",
                        nargs='+',
                        help="output of calculate_methylation_frequency.py",
                        required=True)
    parser.add_argument("-n", "--names",
                        nargs='+',
                        help="names of datasets in --methylation",
                        required=True)
    parser.add_argument("-w", "--window",
                        help="window (region) to which the visualisation has to be restricted",
                        required=True)
    args = parser.parse_args()
    if not len(args.names) == len(args.methylation):
        sys.exit("INPUT ERROR: Expecting the same number of names as datasets!")
    return args


if __name__ == '__main__':
    main()
