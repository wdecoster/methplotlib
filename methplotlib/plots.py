import plotly.graph_objs as go
from methplotlib.annotation import parse_gtf, parse_bed
import sys


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
    annotation = parse_gtf(gtf, window, simplify)
    if annotation:
        for y_pos, transcript in enumerate(annotation):
            line = make_per_gene_line_trace(transcript, window, y_pos)
            exons = [make_per_exon_arrow_trace(transcript, begin, end, y_pos)
                     for begin, end in transcript.exon_tuples
                     if window.begin < begin and window.end > end]
            result.extend([line, *exons])
        return result, y_pos
    else:
        return result, 0


def make_per_gene_line_trace(transcript, window, y_pos):
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


def methylation(meth_data):
    """
    Call function get_data to parse files from nanopolish,
     either the methylation calls (raw) or those from calculate_methylation_frequency
    Return per dataset a list of one (if frequency) or lots of (if raw) plotly traces
    """
    traces = []
    types = []
    names = []
    split = False
    for meth in meth_data:
        if meth.data_type in ['raw', 'phased']:
            traces.append(
                make_per_read_meth_traces(meth.table, phased=meth.data_type == 'phased'))
            split = True
        else:
            traces.append(
                [go.Scatter(x=meth.table.index, y=meth.table["methylated_frequency"],
                            mode='lines',
                            name=meth.name,
                            hoverinfo='name')])
        types.append(meth.data_type)
        names.append(meth.name)
    return DataTraces(traces=traces,
                      types=types,
                      names=names,
                      split=split)


def make_per_read_meth_traces(table, phased=False, max_coverage=100):
    """Make traces for each read"""
    minmax_table = find_min_and_max_pos_per_read(table, phased=phased)
    y_pos_dict = assign_y_height_per_read(minmax_table, phased=phased, max_coverage=max_coverage)
    ratio_cap = min(abs(table["log_lik_ratio"].min()), abs(table["log_lik_ratio"].max()))
    traces = []
    hidden_reads = 0
    for read in table["read_name"].unique():
        strand = table.loc[table["read_name"] == read, "strand"].values[0]
        if phased:
            phase = table.loc[table["read_name"] == read, "HP"].values[0]
        else:
            phase = None
        try:
            traces.append(
                make_per_read_line_trace(read_range=minmax_table.loc[read],
                                         y_pos=y_pos_dict[read],
                                         strand=strand,
                                         phase=phase))
            traces.append(
                make_per_position_likelihood_trace(read_table=table.loc[table["read_name"] == read],
                                                   y_pos=y_pos_dict[read],
                                                   minratio=-ratio_cap,
                                                   maxratio=ratio_cap))
        except KeyError:
            hidden_reads += 1
            continue
    if hidden_reads:
        sys.stderr.write("Warning: hiding {} reads because coverage above {}x.\n".format(
            hidden_reads, max_coverage))
    return traces


def find_min_and_max_pos_per_read(table, phased):
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


def assign_y_height_per_read(df, phased=False, max_coverage=100):
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
    return y_pos


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


def make_per_position_likelihood_trace(read_table, y_pos, minratio, maxratio):
    """Make dots trace per read indicating with the RdBu colorscale from plotly 3.0.0, showing the
    log likelihood of having methylation here"""
    old_RdBu = [[0, 'rgb(5,10,172)'],
                [0.35, 'rgb(106,137,247)'],
                [0.5, 'rgb(190,190,190)'],
                [0.6, 'rgb(220,170,132)'],
                [0.7, 'rgb(230,145,90)'],
                [1, 'rgb(178,10,28)']]
    return go.Scatter(x=read_table['pos'],
                      y=[y_pos] * len(read_table),
                      mode='markers',
                      showlegend=False,
                      text=read_table['log_lik_ratio'],
                      hoverinfo="text",
                      marker=dict(size=4,
                                  color=read_table['log_lik_ratio'],
                                  cmin=minratio,
                                  cmax=maxratio,
                                  colorscale=old_RdBu,
                                  showscale=False))
