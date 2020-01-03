import plotly
import methplotlib.plots as plots
import methplotlib.utils as utils
import methplotlib.qc as qc
from methplotlib.import_methylation import get_data
import logging
import sys


def main():
    args = utils.get_args()
    if args.example:
        utils.print_example()
    utils.init_logs(args)
    windows = utils.make_windows(args.window, fasta=args.fasta)
    for window in windows:
        logging.info("Processing {}".format(window.string))
        meth_data = get_data(args.methylation, args.names, window, args.smooth)
        if args.store:
            import pickle
            pickle.dump(
                obj=meth_data,
                file=open("methplotlib-data-{}.pickle".format(window.string), 'wb'))
        logging.info("Collected methylation data for {} datasets".format(len(meth_data)))
        qc_plots(meth_data, window, qcpath=args.qcfile, outpath=args.outfile)
        logging.info("Created QC plots")
        meth_browser(meth_data=meth_data,
                     window=window,
                     gtf=args.gtf,
                     bed=args.bed,
                     simplify=args.simplify,
                     split=args.split,
                     outfile=args.outfile,
                     dotsize=args.dotsize
                     )
    logging.info("Finished!")


def meth_browser(meth_data, window, gtf=False, bed=False, simplify=False,
                 split=False, outfile=None, dotsize=4):
    """
    meth_Data is a list of Methylation objects from the import_methylation submodule
    annotation is optional and is a gtf or bed file

    if the traces are to be --split per sample (or include raw data as flagged in data.split),
     then show one line per sample and one for the annotation, with methrows = number of datasets
    if no splitting is needed,
     then 4/5 of the browser is used for overlayed samples and one for gtf annotation
    the trace to be used for annotation is thus always num_methrows + 1
    """
    meth_traces = plots.methylation(meth_data, dotsize=dotsize)
    logging.info("Prepared methylation traces.")
    if split or meth_traces.split:
        num_methrows = len(meth_data)
        annot_row = num_methrows + 1
        annot_axis = 'yaxis{}'.format(annot_row)
        fig = create_subplots(num_methrows,
                              split=True,
                              names=meth_traces.names,
                              annotation=bool(bed or gtf))
        for y, (sample_traces, sample_type) in enumerate(meth_traces, start=1):
            logging.info("Adding traces of type {} at height {}".format(sample_type, y))
            for meth_trace in sample_traces:
                fig.append_trace(trace=meth_trace, row=y, col=1)
            if sample_type == 'nanopolish_freq':
                fig["layout"]["yaxis{}".format(y)].update(title="Modified <br> frequency")
                fig["layout"].update(showlegend=False)
                fig["layout"]["xaxis"].update(tickformat='g',
                                              separatethousands=True,
                                              range=[window.begin, window.end])
                fig["layout"].update(legend=dict(orientation='h'))
            elif sample_type in ['nanopolish_call', 'nanopolish_phased']:
                fig["layout"]["yaxis{}".format(y)].update(title="Reads")
                fig["layout"].update(showlegend=False)
                fig["layout"]["xaxis"].update(tickformat='g',
                                              separatethousands=True,
                                              range=[window.begin, window.end])
            elif sample_type == 'nanocompore':
                fig["layout"]["yaxis{}".format(y)].update(title="-log10(pval)")
                fig["layout"]["xaxis"].update(tickformat='g',
                                              range=[window.begin, window.end])
                fig["layout"].update(legend=dict(orientation='h'))
            elif sample_type == 'ont-cram':
                fig["layout"]["yaxis{}".format(y)].update(title="Reads")
                fig["layout"]["xaxis"].update(tickformat='g',
                                              separatethousands=True,
                                              range=[window.begin, window.end])
            else:
                sys.exit("ERROR: unrecognized data type {}".format(sample_type))
    else:
        num_methrows = 4
        annot_row = 5
        annot_axis = 'yaxis2'
        fig = create_subplots(num_methrows, split=False, annotation=bool(bed or gtf))
        for meth_trace in meth_traces.traces:
            fig.append_trace(trace=meth_trace[0], row=1, col=1)
        fig["layout"].update(legend=dict(orientation='h'))
        fig["layout"]['yaxis'].update(title="Modified <br> frequency")
        fig["layout"]["xaxis"].update(tickformat='g',
                                      separatethousands=True,
                                      range=[window.begin, window.end])
    logging.info("Prepared modification plots.")
    if bed:
        for annot_trace in plots.bed_annotation(bed, window):
            fig.append_trace(trace=annot_trace, row=annot_row, col=1)
        y_max = -2
    if gtf:
        annotation_traces, y_max = plots.gtf_annotation(gtf, window, simplify)
        for annot_trace in annotation_traces:
            fig.append_trace(trace=annot_trace, row=annot_row, col=1)
    if bed or gtf:
        fig["layout"][annot_axis].update(range=[-2, y_max + 1],
                                         showgrid=False,
                                         zeroline=False,
                                         showline=False,
                                         ticks='',
                                         showticklabels=False)
        logging.info("Prepared annotation plots.")
    fig["layout"].update(barmode='overlay',
                         title="Nucleotide modifications",
                         hovermode='closest',
                         plot_bgcolor='rgba(0,0,0,0)')
    utils.create_browser_output(fig, outfile, window)


def create_subplots(num_methrows, split, names=None, annotation=True):
    '''
    Prepare the panels (rows * 1 column) for the subplots.
    If splitting: one row for each dataset, taking 90%/len(datasets) for heights
    If not: one row spanning 4 rows and taking 90% of the heights
    if annotation is True (bed or gtf) then add a row with height 10%
    '''
    if split:
        return plotly.subplots.make_subplots(
            rows=num_methrows + annotation,
            cols=1,
            shared_xaxes=True,
            specs=[[{}] for i in range(num_methrows + annotation)],
            print_grid=False,
            subplot_titles=names,
            vertical_spacing=0.1,
            row_heights=[0.9 / num_methrows] * num_methrows + [0.1] * annotation

        )
    else:
        return plotly.subplots.make_subplots(
            rows=num_methrows + annotation,
            cols=1,
            shared_xaxes=True,
            specs=[[{'rowspan': num_methrows}], [None], [None], [None]] + [[{}]] * annotation,
            print_grid=False,
            vertical_spacing=0.1,
            row_heights=[0.9, 0, 0, 0] + [0.1] * annotation
        )


def qc_plots(meth_data, window, qcpath=None, outpath=None):

    if qcpath is None and outpath is None:
        outfile = "qc_report_{}.html".format(window.string)
    elif qcpath is None:
        from pathlib import Path, PosixPath
        p = Path(outpath)
        Path.mkdir(p.parent, exist_ok=True, parents=True)
        outfile = str(p.parent / PosixPath("qc_" + p.stem + ".html"))
    else:
        from pathlib import Path
        p = Path(qcpath)
        Path.mkdir(p.parent, exist_ok=True, parents=True)
        outfile = qcpath

    with open(outfile, 'w') as qc_report:
        qc_report.write(qc.num_sites_bar(meth_data))
        if len([m for m in meth_data if m.data_type == "nanopolish_freq"]) > 0:
            data = [m.table.rename({"methylated_frequency": m.name}, axis='columns')
                    for m in meth_data if m.data_type == "nanopolish_freq"]
            full = data[0].join(data[1:]).dropna(how="any", axis="index")
            qc_report.write(qc.modified_fraction_histogram(full))
        if len([m for m in meth_data if m.data_type == "nanopolish_freq"]) > 2:
            qc_report.write(qc.pairwise_correlation_plot(full))
            qc_report.write(qc.pca(full))
            qc_report.write(qc.global_box(data))
        if len([m for m in meth_data
                if m.data_type in ["nanopolish_call", "nanopolish_phased"]]) > 2:
            pass


if __name__ == '__main__':
    main()
