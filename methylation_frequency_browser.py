import pandas as pd
import numpy as np
import sys
import plotly
from plotly import tools
import plotly.graph_objs as go
from argparse import ArgumentParser
from gtfparse import read_gtf
import itertools


class Transcript(object):
    def __init__(self, transcript_id, name, exon_tuples, strand):
        self.transcript_id = transcript_id
        self.name = name
        self.exon_tuples = list(exon_tuples)
        self.intron_tuples = self.exon_tuples_to_intron_tuples()
        self.strand = strand
        self.begin = min(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.end = max(list(itertools.chain.from_iterable(self.exon_tuples)))

    def exon_tuples_to_intron_tuples(self):
        """Convert a list of exon tuples to intron tuples
        Flattens the list, drops the first and last coordinate and creates new tuples
        """
        introns_coords = iter(list(itertools.chain.from_iterable(self.exon_tuples))[1:-1])
        return list(zip(introns_coords, introns_coords))


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
    """
    Import methylation frequency from all files in args.methylation
    within the window args.window
    """
    return [read_meth_freq(f, args.window) for f in args.methylation]


def parse_gtf(gtff, window):
    """
    Parse the gtff using read_gtf and select the relevant region
    as determined by the window

    """
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
        sys.stderr.write("Found {} overlapping transcripts.\n".format(len(result)))
        return result


def plotly_annotation(annotation, window):
    """
    Return a plotly trace for the annotation
    with a line for the entire gene and thicker bars for exons
    """
    chromosome, begin, end = parse_region(window)
    result = []
    for y_pos, transcript in enumerate(annotation):
        line = go.Scatter(x=[max(transcript.begin, begin), min(transcript.end, end)],
                          y=[y_pos, y_pos],
                          mode='lines',
                          line=dict(width=2),
                          name=transcript.transcript_id,
                          text=transcript.name,
                          hoverinfo='text',
                          showlegend=False)
        exons = [go.Scatter(x=[begin, end],
                            y=[y_pos, y_pos],
                            mode='lines',
                            line=dict(width=8),
                            name=transcript.transcript_id,
                            text=transcript.name,
                            hoverinfo='text',
                            showlegend=False)
                 for begin, end in transcript.exon_tuples]
        result.extend([line, *exons])
    return result, y_pos


def plotly_methylation(meth, name):
    """
    Return a plotly trace for the methylation frequency of this dataset
    """
    return go.Scatter(x=meth.index, y=meth["methylated_frequency"],
                      mode='lines+markers',
                      name=name,
                      hoverinfo='name')


def meth_browser(methlist, names, window, annotation=False):
    """
    methlist is a list of pandas dataframes containing 'chromosome', 'pos', 'methylated_frequency'
    names should have the same length as methlist and contain identifiers for the datasets
    annotation is optional and is a gtf processed by parse_gtf()
    """
    fig = tools.make_subplots(rows=5,
                              cols=1,
                              shared_xaxes=True,
                              specs=[
                                  [{'rowspan': 4}],
                                  [None],
                                  [None],
                                  [None],
                                  [{}],
                              ],
                              print_grid=False
                              )
    for meth_trace in [plotly_methylation(a, n) for a, n in zip(methlist, names)]:
        fig.append_trace(trace=meth_trace,
                         row=1,
                         col=1)
    if annotation:
        annotation_traces, y_max = plotly_annotation(annotation, window)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace,
                             row=5,
                             col=1)
    fig["layout"].update(barmode='overlay',
                         title="Methylation frequency",
                         hovermode='closest')
    fig["layout"]["yaxis5"].update(range=[0, y_max],
                                   showgrid=False,
                                   zeroline=False,
                                   showline=False,
                                   ticks='',
                                   showticklabels=False)
    chromosome, begin, end = parse_region(window)
    fig["layout"]["xaxis5"].update(range=[begin, end])
    html = plotly.offline.plot(fig,
                               output_type="div",
                               show_link=False)
    with open("methylation_browser.html", 'w') as output:
        output.write(html)


def main():
    args = get_args()
    methlist = get_data(args)
    annotation = parse_gtf(args.gtf, args.window)
    meth_browser(methlist, names=args.names, window=args.window, annotation=annotation)


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
