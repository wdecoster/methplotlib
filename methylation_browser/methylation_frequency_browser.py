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
        self.strand = strand
        self.begin = min(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.end = max(list(itertools.chain.from_iterable(self.exon_tuples)))
        self.color = ""


class Region(object):
    def __init__(self, region):
        self.chromosome, interval = region.split(':')
        self.begin, self.end = [int(i) for i in interval.split('-')]


def read_meth_freq(filename, window):
    table = pd.read_csv(filename, sep="\t")
    table = table.loc[(table.chromosome == window.chromosome) &
                      (table.start > window.begin) &
                      (table.end < window.end)]
    table["pos"] = np.floor(table[['start', 'end']].mean(axis=1))
    table = table.drop(columns=['start', 'end', 'num_cpgs_in_group',
                                'called_sites', 'called_sites_methylated', 'group_sequence']) \
        .sort_values('pos') \
        .groupby('pos') \
        .mean() \
        .rolling(window=5, center=True) \
        .mean()
    return table


def get_data(methylation_files, window):
    """
    Import methylation frequency from all files in args.methylation
    within the window args.window
    """
    return [read_meth_freq(f, window) for f in methylation_files]


def parse_gtf(gtff, window):
    """
    Parse the gtff using read_gtf and select the relevant region
    as determined by the window

    """
    if not gtff:
        return False
    else:
        gtf = read_gtf(gtff)
        columns = ["start", "end", "strand", "transcript_id", "gene_name"]
        gtf_f = gtf.loc[(gtf["feature"] == "exon") & (gtf["seqname"] == window.chromosome), columns]
        transcript_slice = (gtf_f.groupby("transcript_id")["start"].max() > window.begin) & (
            gtf_f.groupby("transcript_id")["end"].min() < window.end)
        transcripts = transcript_slice[transcript_slice].index
        region = gtf_f.loc[gtf_f["transcript_id"].isin(transcripts)]
        result = []
        for t in transcripts:
            tr = region.loc[region["transcript_id"] == t]
            result.append(
                Transcript(transcript_id=t,
                           name=tr["gene_name"].tolist()[0],
                           exon_tuples=tr.loc[:, ["start", "end"]]
                           .sort_values("start")
                           .itertuples(index=False, name=None),
                           strand=tr["strand"].tolist()[0])
            )
        sys.stderr.write("Found {} overlapping transcripts.\n".format(len(result)))
        genes = set([t.name for t in result])
        colors = plotly.colors.DEFAULT_PLOTLY_COLORS * 100
        colordict = {g: c for g, c in zip(genes, colors)}
        for t in result:
            t.color = colordict[t.name]
        return result


def plotly_annotation(annotation, window):
    """
    Return a plotly trace for the annotation
    with a line for the entire gene and thicker bars for exons
    """
    result = []
    for y_pos, transcript in enumerate(annotation):
        line = go.Scatter(x=[max(transcript.begin, window.begin), min(transcript.end, window.end)],
                          y=[y_pos, y_pos],
                          mode='lines',
                          line=dict(width=2, color=transcript.color),
                          name=transcript.transcript_id,
                          text=transcript.name,
                          hoverinfo='text',
                          showlegend=False)
        if transcript.strand == "+":
            marker = "triangle-right"
        else:
            marker = "triangle-left"
        exons = [go.Scatter(x=[begin, end],
                            y=[y_pos, y_pos],
                            mode='lines+markers',
                            line=dict(width=8, color=transcript.color),
                            name=transcript.transcript_id,
                            text=transcript.name,
                            hoverinfo='text',
                            showlegend=False,
                            marker=dict(symbol=marker,
                                        size=8))
                 for begin, end in transcript.exon_tuples
                 if window.begin < begin and window.end > end]
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
        fig["layout"]["yaxis5"].update(range=[0, y_max],
                                       showgrid=False,
                                       zeroline=False,
                                       showline=False,
                                       ticks='',
                                       showticklabels=False)
        fig["layout"]["xaxis5"].update(range=[window.begin, window.end])
    fig["layout"].update(barmode='overlay',
                         title="Methylation frequency",
                         hovermode='closest')
    html = plotly.offline.plot(fig,
                               output_type="div",
                               show_link=False)
    with open("methylation_browser.html", 'w') as output:
        output.write(html)


def main():
    args = get_args()
    window = Region(args.window)
    methlist = get_data(args.methylation, window)
    annotation = parse_gtf(args.gtf, window)
    meth_browser(methlist, names=args.names, window=window, annotation=annotation)


def get_args():
    parser = ArgumentParser(description="plotting methylation frequency")
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
