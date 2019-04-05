import plotly
from plotly import tools
import methplotlib.plots as plots
import methplotlib.utils as utils


def main():
    args = utils.get_args()
    meth_browser(methlist=args.methylation,
                 names=args.names,
                 window=utils.Region(args.window),
                 gtf=args.gtf,
                 smoothen=args.smooth)


def meth_browser(methlist, names, window, gtf=False, smoothen=5):
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
    for meth_trace in plots.methylation(methlist, names, window, smoothen):
        fig.append_trace(trace=meth_trace, row=1, col=1)
    if gtf:
        annotation_traces, y_max = plots.annotation(gtf, window)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace, row=5, col=1)
        fig["layout"]["yaxis2"].update(range=[0, y_max],
                                       showgrid=False,
                                       zeroline=False,
                                       showline=False,
                                       ticks='',
                                       showticklabels=False)
        fig["layout"]["xaxis"].update(range=[window.begin, window.end])
    fig["layout"].update(barmode='overlay',
                         title="Methylation frequency",
                         hovermode='closest')
    with open("methylation_browser.html", 'w') as output:
        output.write(plotly.offline.plot(fig,
                                         output_type="div",
                                         show_link=False))


if __name__ == '__main__':
    main()
