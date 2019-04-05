import plotly.graph_objs as go
from methplotlib.annotation import parse_gtf
from methplotlib.import_methylation import get_data


def annotation(gtf, window):
    """
    Return a plotly trace for the annotation
    with a line for the entire gene and triangles for exons,
    indicating direction of transcription
    """
    result = []
    annotation = parse_gtf(gtf, window)
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


def methylation(methylation_files, names, window, smoothen):
    """
    Call function get_data to parse files from calculate_methylation_frequency
    Return a plotly trace for the methylation frequency of this dataset
    """
    meth = get_data(methylation_files, window, smoothen)
    return [go.Scatter(x=meth.index, y=meth["methylated_frequency"],
                       mode='lines',
                       name=name,
                       hoverinfo='name')
            for meth, name in zip(meth, names)]
