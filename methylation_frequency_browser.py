import pandas as pd
import numpy as np
import sys
import plotly
import plotly.graph_objs as go
from argparse import ArgumentParser
from gtfparse import read_gtf


class Transcript(object):
    def __init__(self, identifier, name, exontuples, strand):
        self.identifier = identifier
        self.name = name
        self.exontuples = list(exontuples)
        self.strand = strand


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


def parse_gtf(gtff, window):
    if not gtff:
        return False
    else:
        chromosome, begin, end = parse_region(window)
        gtf = read_gtf(gtff)
        columns = ["start", "end", "strand", "transcript_id", "transcript_name"]
        gtf_f = gtf.loc[(gtf["feature"] == "exon") & (gtf["seqname"] == chromosome), columns]
        transcript_slice = (gtf_f.groupby("transcript_id")["start"].max() > begin) & (
            gtf_f.groupby("transcript_id")["end"].min() < end)
        transcripts = transcript_slice[transcript_slice].index
        region = gtf_f.loc[gtf_f["transcript_id"].isin(transcripts)]
        result = []
        for t in transcripts:
            tr = region.loc[region["transcript_id"] == t]
            result.append(
                Transcript(transcript_id=t,
                           name=tr["transcript_name"].tolist()[0],
                           exon_tuples=tr.loc[:, ["start", "end"]]
                           .sort_values("start")
                           .itertuples(index=False, name=None),
                           strand=tr["strand"].tolist()[0])
            )


def meth_browser(methlist, names, annotation=False):
    """
    methlist is a list of pandas dataframes containing 'chromosome', 'pos', 'methylated_frequency'
    names should have the same length as methlist and contain identifiers for the datasets
    annotation is optional and is a gtf processed by parse_gtf()
    """
    data = [go.Scatter(x=a.index, y=a["methylated_frequency"],
                       mode='lines+markers',
                       name=n) for a, n in zip(methlist, names)]
    if annotation:
        data = data + plotly_gtf(annotation)
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
    annotation = parse_gtf(args.gtf, args.window)
    meth_browser(methlist, names=args.names, annotation=annotation)


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
    parser.add_argument("-g", "--gtf",
                        help="add annotation based on a gtf file matching to your reference genome")
    args = parser.parse_args()
    if not len(args.names) == len(args.methylation):
        sys.exit("INPUT ERROR: Expecting the same number of names as datasets!")
    return args


if __name__ == '__main__':
    main()
