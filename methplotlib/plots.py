import plotly.graph_objs as go
from methplotlib.annotation import parse_gtf


class DataTraces(object):
    def __init__(self, traces, split):
        self.traces = traces
        self.split = split


def annotation(gtf, window, simplify=False):
    """
    Return a plotly trace for the annotation
    with a line for the entire gene and triangles for exons,
    indicating direction of transcription
    """
    result = []
    annotation = parse_gtf(gtf, window, simplify)
    if annotation:
        for y_pos, transcript in enumerate(annotation):
            line = gene_line_trace(transcript, window, y_pos)
            exons = [exon_arrow_trace(transcript, begin, end, y_pos)
                     for begin, end in transcript.exon_tuples
                     if window.begin < begin and window.end > end]
            result.extend([line, *exons])
        return result, y_pos
    else:
        return result, 0


def gene_line_trace(transcript, window, y_pos):
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


def exon_arrow_trace(transcript, begin, end, y_pos):
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


def methylation(meth_data):
    """
    Call function get_data to parse files from nanopolish,
     either the methylation calls (raw) or those from calculate_methylation_frequency
    Return per dataset a list of one (if frequency) or lots of (if raw) plotly traces
    """
    traces = []
    split = False
    for meth in meth_data:
        if meth.data_type in ['raw', 'phased']:
            traces.append(per_read_traces(meth.table, phased=meth.data_type == 'phased'))
            split = True
        else:
            traces.append([go.Scatter(x=meth.table.index, y=meth.table["methylated_frequency"],
                                      mode='lines',
                                      name=meth.name,
                                      hoverinfo='name')])
    return DataTraces(traces=traces, split=split)


def per_read_traces(table, phased=False):
    """Make traces for each read"""
    minmax_table = min_and_max_pos(table, phased=phased)
    y_pos_dict = assign_y_pos(minmax_table, phased=phased)
    ratio_cap = min(abs(table["log_lik_ratio"].min()), abs(table["log_lik_ratio"].max()))
    traces = []
    for read in table["read_name"].unique():
        strand = table.loc[table["read_name"] == read, "strand"].values[0]
        if phased:
            phase = table.loc[table["read_name"] == read, "HP"].values[0]
        else:
            phase = None
        traces.append(read_line_trace(read_range=minmax_table.loc[read],
                                      y_pos=y_pos_dict[read],
                                      strand=strand,
                                      phase=phase))
        traces.append(position_likelihood_trace(read_table=table.loc[table["read_name"] == read],
                                                y_pos=y_pos_dict[read],
                                                minratio=-ratio_cap,
                                                maxratio=ratio_cap))
    return traces


def min_and_max_pos(table, phased):
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


def assign_y_pos(df, phased=False):
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
    heights = [[] for i in range(100)]
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


def read_line_trace(read_range, y_pos, strand, phase=None):
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


def position_likelihood_trace(read_table, y_pos, minratio, maxratio):
    """Make dots trace per read indicating (with RdBu) colorscale the
    log likelihood of having methylation here"""
    return go.Scatter(x=read_table['pos'],
                      y=[y_pos] * len(read_table),
                      mode='markers',
                      showlegend=False,
                      marker=dict(size=4,
                                  color=read_table['log_lik_ratio'],
                                  cmin=minratio,
                                  cmax=maxratio,
                                  colorscale='RdBu',
                                  showscale=True))


def splom(meth_data):
    data = [m.table.rename({"methylated_frequency": m.name}, axis='columns')
            for m in meth_data if m.data_type == "frequency"]
    labels = [m.name for m in meth_data if m.data_type == "frequency"]
    full = data[0].join(data[1:])
    trace = go.Splom(dimensions=[dict(label=l, values=full[l]) for l in labels],
                     marker=dict(size=4,
                                 line=dict(width=0.5,
                                           color='rgb(230,230,230)')),
                     diagonal=dict(visible=False)
                     )

    layout = go.Layout(
        title='Correlation of methylation frequency',
        dragmode='select',
        width=1200,
        height=1200,
        autosize=False,
        hovermode=False,
        plot_bgcolor='rgba(240,240,240, 0.95)')

    for i in len(data):
        layout["xaxis{}".format(i)] = dict(
            showline=True, zeroline=False, gridcolor='#fff', ticklen=4)
        layout["yaxis{}".format(i)] = dict(
            showline=True, zeroline=False, gridcolor='#fff', ticklen=4)

    return dict(data=[trace], layout=layout)
