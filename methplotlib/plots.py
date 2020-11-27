import plotly.graph_objs as go
from methplotlib.annotation import parse_annotation, parse_bed
import sys
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import numpy as np


class DataTraces(object):
    def __init__(self, traces, types, names, split):
        self.traces = traces
        self.types = types
        self.names = names
        self.split = split
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.index == len(self.traces):
            raise StopIteration
        else:
            self.index += 1
            return self.traces[self.index - 1], self.types[self.index - 1]


def gtf_annotation(gtf, window, simplify=False):
    """
    Return a plotly trace for the annotation
    with a line for the entire gene and triangles for exons,
    indicating direction of transcription
    """
    result = []
    annotation = parse_annotation(gtf, window, simplify)
    if annotation:
        for y_pos, transcript in enumerate(annotation):
            line = make_per_gene_annot_line_trace(transcript, window, y_pos)
            exons = [make_per_exon_arrow_trace(transcript, begin, end, y_pos)
                     for begin, end in transcript.exon_tuples
                     if window.begin < begin and window.end > end]
            result.extend([line, *exons])
        return result, y_pos
    else:
        return result, 0


def make_per_gene_annot_line_trace(transcript, window, y_pos):
    """Generate a line trace for the gene

    Trace can get limited by the window sizes
    """
    return go.Scatter(x=[max(transcript.begin, window.begin),
                         min(transcript.end, window.end)],
                      y=[y_pos, y_pos],
                      mode='lines',
                      line=dict(width=2, color=transcript.color),
                      name=transcript.transcript,
                      text=transcript.gene,
                      hoverinfo='text',
                      showlegend=False)


def make_per_exon_arrow_trace(transcript, begin, end, y_pos):
    """Generate a line+marker trace for the exon

    The shape is an arrow, as defined by the strand in transcript.marker
    """
    return go.Scatter(x=[begin, end],
                      y=[y_pos, y_pos],
                      mode='lines+markers',
                      line=dict(width=8, color=transcript.color),
                      name=transcript.transcript,
                      text=transcript.gene,
                      hoverinfo='text',
                      showlegend=False,
                      marker=dict(symbol=transcript.marker,
                                  size=8))


def bed_annotation(bed, window):
    return [go.Scatter(x=[begin, end],
                       y=[-2, -2],
                       mode='lines',
                       line=dict(width=16, color='grey'),
                       text=name,
                       hoverinfo='text',
                       showlegend=False)
            for (begin, end, name) in parse_bed(bed, window)]


def methylation(meth_data, dotsize=4, binary=False):
    """
    Plot methylation traces from various data types
    """
    traces = []
    types = []
    names = []
    split = False
    for meth in meth_data:
        if meth.data_type in ['nanopolish_call', 'nanopolish_phased']:
            traces.append(
                make_per_read_meth_traces_llr(table=meth.table,
                                              phased=meth.data_type == 'nanopolish_phased',
                                              dotsize=dotsize,
                                              binary=binary)
            )
            split = True
        elif meth.data_type == 'nanocompore':
            traces.append(plot_nanocompore(meth.table, dotsize=dotsize))
            split = True
        elif meth.data_type == 'nanopolish_freq':
            traces.append(
                [go.Scatter(x=meth.table.index, y=meth.table["methylated_frequency"],
                            mode='lines',
                            name=meth.name,
                            hoverinfo='name')])
        elif meth.data_type == 'ont-cram':
            traces.append(
                make_per_read_meth_traces_phred(table=meth.table,
                                                dotsize=dotsize)
            )
            split = True
        elif meth.data_type == 'bedgraph':
            starts = meth.table["Start"]
            ends = meth.table["End"]
            vals = meth.table["Value"]
            traces.append(
                [go.Scatter(x=[val for pair in zip(starts, ends) for val in pair],
                            y=[val for pair in zip(vals, vals) for val in pair],
                            mode='lines',
                            name=meth.name,
                            hoverinfo='name')])
        else:
            sys.exit(f"ERROR: unrecognized data type {meth.data_type}")
        types.append(meth.data_type)
        names.append(meth.name)
    return DataTraces(traces=traces,
                      types=types,
                      names=names,
                      split=split)


def make_per_read_meth_traces_phred(table, max_cov=100, dotsize=4):
    """Make traces for each read"""
    minmax_table = find_min_and_max_pos_per_read(table)
    df_heights = assign_y_height_per_read(minmax_table, max_coverage=max_cov)
    table = table.join(df_heights, on="read_name")
    traces = []
    hidden = 0
    for read in table["read_name"].unique():
        strand = table.loc[table["read_name"] == read, "strand"].values[0]
        try:
            traces.append(
                make_per_read_line_trace(read_range=minmax_table.loc[read],
                                         y_pos=df_heights.loc[read, 'height'],
                                         strand=strand)
            )
        except KeyError:
            hidden += 1
            continue
    if hidden:
        sys.stderr.write(f"Warning: hiding {hidden} reads because coverage above {max_cov}x.\n")
    traces.append(
        make_per_position_phred_scatter(read_table=table, dotsize=dotsize)
    )
    return traces


def make_per_read_meth_traces_llr(table, phased=False, max_cov=100, dotsize=4, binary=False):
    """Make traces for each read"""
    minmax_table = find_min_and_max_pos_per_read(table, phased=phased)
    df_heights = assign_y_height_per_read(minmax_table, phased=phased, max_coverage=max_cov)
    table = table.join(df_heights, on="read_name")
    if binary:
        table.loc[:, "llr_scaled"] = binarize_log_likelihood_ratio(table["log_lik_ratio"].copy())
    else:
        table.loc[:, "llr_scaled"] = rescale_log_likelihood_ratio(table["log_lik_ratio"].copy())
    traces = []
    hidden = 0
    for read in table["read_name"].unique():
        strand = table.loc[table["read_name"] == read, "strand"].values[0]
        if phased:
            phase = table.loc[table["read_name"] == read, "HP"].values[0]
        else:
            phase = None
        try:
            traces.append(
                make_per_read_line_trace(read_range=minmax_table.loc[read],
                                         y_pos=df_heights.loc[read, 'height'],
                                         strand=strand,
                                         phase=phase)
            )
        except KeyError:
            hidden += 1
            continue
    if hidden:
        sys.stderr.write(f"Warning: hiding {hidden} reads because coverage above {max_cov}x.\n")
    traces.append(
        make_per_position_likelihood_scatter(read_table=table, dotsize=dotsize)
    )
    return traces


def find_min_and_max_pos_per_read(table, phased=False):
    """Return a table with for every read the minimum and maximum position"""
    mm_table = table.loc[:, ["read_name", "pos"]] \
        .groupby('read_name') \
        .min() \
        .join(table.loc[:, ["read_name", "pos"]]
              .groupby('read_name')
              .max(),
              lsuffix="min",
              rsuffix="max")
    if phased:
        return mm_table.join(table.loc[:, ["read_name", "HP"]]
                             .drop_duplicates(subset='read_name')
                             .set_index('read_name'))
    else:
        return mm_table


def assign_y_height_per_read(df, phased=False, max_coverage=1000):
    """Assign height of the read in the per read traces

    Gets a dataframe of read_name, posmin and posmax.
    Sorting by position, and optionally by phase block.
    Determines optimal height (y coordinate) for this read

    Returns a dictionary mapping read_name to y_coord
    """
    if phased:
        dfs = df.sort_values(by=['HP', 'posmin', 'posmax'],
                             ascending=[True, True, False])
    else:
        dfs = df.sort_values(by=['posmin', 'posmax'],
                             ascending=[True, False])
    heights = [[] for i in range(max_coverage)]
    y_pos = dict()
    for read in dfs.itertuples():
        for y, layer in enumerate(heights, start=1):
            if len(layer) == 0:
                layer.append(read.posmax)
                y_pos[read.Index] = y
                break
            if read.posmin > layer[-1]:
                layer.append(read.posmax)
                y_pos[read.Index] = y
                break
    return pd.DataFrame({'read': list(y_pos.keys()),
                         'height': list(y_pos.values())}) \
        .set_index('read')


def rescale_log_likelihood_ratio(llr):
    """
    Rescale log likelihood ratios

    positive ratios between 0 and 1
    negative ratios between -1 and 0
    """
    scaler = MinMaxScaler(feature_range=(0, 1))
    llr[llr > 0] = scaler.fit_transform(llr[llr > 0].values.reshape(-1, 1)).tolist()
    scaler = MinMaxScaler(feature_range=(-1, 0))
    llr[llr < 0] = scaler.fit_transform(llr[llr < 0].values.reshape(-1, 1)).tolist()
    return llr


def binarize_log_likelihood_ratio(llr, cutoff=2):
    '''
    Converts the log likelihood to sorta binary representation,
    with 0 being insufficient evidence (based on cutoff)
    and +1 being probably methylated
    and -1 being probably unmethylated
    '''
    llr[(-cutoff < llr) & (llr < cutoff)] = 0
    llr[llr > cutoff] = 1
    llr[llr < -cutoff] = -1
    return llr


def make_per_read_line_trace(read_range, y_pos, strand, phase=None):
    """Make a grey line trace for a single read,
    with black arrow symbols on the edges indicating strand"""
    symbol = "triangle-right" if strand == "+" else "triangle-left"
    if phase:
        if phase == 1:
            color = 'lightgreen'
        elif phase == 2:
            color = 'yellow'
        else:  # phase is np.nan
            color = 'black'
    else:
        color = 'black'
    return go.Scatter(x=[read_range['posmin'], read_range['posmax']],
                      y=[y_pos, y_pos],
                      mode='lines+markers',
                      line=dict(width=1, color='lightgrey'),
                      showlegend=False,
                      marker=dict(symbol=symbol,
                                  size=8,
                                  color=color,
                                  line=dict(width=0.5,
                                            color='black')))


def make_per_position_likelihood_scatter(read_table, maxval=0.75, dotsize=4):
    """Make scatter plot per CpG per read
    with the RdBu colorscale from plotly 3.0.0, showing the
    scaled log likelihood of having methylation here"""
    old_RdBu = [[0, 'rgb(5,10,172)'],
                [0.35, 'rgb(106,137,247)'],
                [0.5, 'rgb(190,190,190)'],
                [0.6, 'rgb(220,170,132)'],
                [0.7, 'rgb(230,145,90)'],
                [1, 'rgb(178,10,28)']]
    return go.Scatter(x=read_table['pos'],
                      y=read_table['height'],
                      mode='markers',
                      showlegend=False,
                      text=read_table['log_lik_ratio'],
                      hoverinfo="text",
                      marker=dict(size=dotsize,
                                  color=read_table['llr_scaled'],
                                  cmin=-maxval,
                                  cmax=maxval,
                                  colorscale=old_RdBu,
                                  showscale=True,
                                  colorbar=dict(title="Modification likelihood",
                                                titleside="right",
                                                tickvals=[-maxval, 0, maxval],
                                                ticktext=["Likely <br> unmodified",
                                                          "0", "Likely <br> modified"],
                                                ticks="outside"))
                      )


def make_per_position_phred_scatter(read_table, dotsize=4):
    """Make scatter plot per CpG per read"""
    return go.Scatter(x=read_table['pos'],
                      y=read_table['height'],
                      mode='markers',
                      showlegend=False,
                      text=read_table['mod'],
                      hoverinfo="text",
                      marker=dict(size=dotsize,
                                  color=read_table['quality'],
                                  colorscale='Reds',
                                  colorbar=dict(title="Modification probability",
                                                titleside="right",
                                                tickvals=[read_table['quality'].min(),
                                                          read_table['quality'].max()],
                                                ticktext=["Likely <br> unmodified",
                                                          "Likely <br> modified"],
                                                ticks="outside")
                                  ))


def plot_nanocompore(table, dotsize=4):
    return [go.Scatter(x=table['pos'],
                       y=-table[pval].apply(np.log10),
                       mode='markers+lines',
                       name=pval,
                       marker=dict(size=dotsize))
            for pval in table.columns if not pval == 'pos']
