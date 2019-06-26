import plotly.graph_objs as go
from methplotlib.annotation import parse_gtf
from methplotlib.import_methylation import get_data


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


def methylation(methylation_files, names, window, smoothen):
    """
    Call function get_data to parse files from nanopolish,
     either the methylation calls (raw) or those from calculate_methylation_frequency
    Return a plotly trace for the methylation visualisation of this dataset
    """
    traces = []
    split = False
    for meth, name in zip(get_data(methylation_files, window, smoothen), names):
        if meth.data_type == 'raw':
            traces.append()
            split = True
        else:
            traces.append(go.Scatter(x=meth.index, y=meth["methylated_frequency"],
                                     mode='lines',
                                     name=name,
                                     hoverinfo='name'))
    return DataTraces(traces=traces, split=split)
