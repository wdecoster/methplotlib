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
        logging.info(f"Processing {window.string}")
        meth_data = get_data(args.methylation, args.names, window, args.smooth)
        if args.store:
            import pickle
            pickle.dump(
                obj=meth_data,
                file=open(f"methplotlib-data-{window.string}.pickle", 'wb'))
        logging.info(f"Collected methylation data for {len(meth_data)} datasets")
        qc.qc_plots(meth_data, window, qcpath=args.qcfile, outpath=args.outfile)
        logging.info("Created QC plots")
        meth_browser(meth_data=meth_data,
                     window=window,
                     gtf=args.gtf,
                     bed=args.bed,
                     simplify=args.simplify,
                     split=args.split,
                     outfile=args.outfile,
                     dotsize=args.dotsize,
                     static=args.static,
                     binary=args.binary,
                     )
    logging.info("Finished!")


def meth_browser(meth_data, window, gtf=False, bed=False, simplify=False,
                 split=False, outfile=None, dotsize=4, static=False, binary=False):
    """
    meth_Data is a list of Methylation objects from the import_methylation submodule
    annotation is optional and is a gtf or bed file

    if the traces are to be --split per sample (or include raw data as flagged in data.split),
     then show one line per sample and one for the annotation, with methrows = number of datasets
    if no splitting is needed,
     then 4/5 of the browser is used for overlayed samples and one for gtf annotation
    the trace to be used for annotation is thus always num_methrows + 1
    """
    meth_traces = plots.methylation(meth_data, dotsize=dotsize, binary=binary)
    logging.info("Prepared methylation traces.")
    if split or meth_traces.split:
        num_methrows = len(meth_data)
        logging.info(f'Making browser in split mode, with {num_methrows} modification rows.')
        annot_row = num_methrows + 1
        annot_axis = f'yaxis{annot_row}'
        fig = utils.create_subplots(num_methrows,
                                    split=True,
                                    names=meth_traces.names,
                                    annotation=bool(bed or gtf))
        for y, (sample_traces, sample_type) in enumerate(meth_traces, start=1):
            logging.info(f"Adding traces of type {sample_type} at height {y}")
            for meth_trace in sample_traces:
                fig.append_trace(trace=meth_trace, row=y, col=1)
            if sample_type == 'nanopolish_freq':
                fig["layout"][f"yaxis{y}"].update(title="Modified <br> frequency")
                fig["layout"].update(showlegend=False)
                fig["layout"].update(legend=dict(orientation='h'))
            elif sample_type in ['nanopolish_call', 'nanopolish_phased']:
                fig["layout"][f"yaxis{y}"].update(title="Reads")
                fig["layout"].update(showlegend=False)
            elif sample_type == 'nanocompore':
                fig["layout"][f"yaxis{y}"].update(title="-log10(pval)")
                fig["layout"].update(legend=dict(orientation='h'))
            elif sample_type == 'ont-cram':
                fig["layout"][f"yaxis{y}"].update(title="Reads")
            elif sample_type == 'bedgraph':
                fig["layout"][f"yaxis{y}"].update(title="Value")
            else:
                sys.exit(f"ERROR: unrecognized data type {sample_type}")
    else:
        logging.info('Making browser in overlaying mode.')
        num_methrows = 4
        annot_row = 5
        annot_axis = 'yaxis2'
        fig = utils.create_subplots(num_methrows, split=False, annotation=bool(bed or gtf))
        for meth_trace in meth_traces.traces:
            for trace in meth_trace:
                fig.append_trace(trace=trace, row=1, col=1)
        fig["layout"].update(legend=dict(orientation='h'))
        if meth_traces.types[0] == 'nanopolish_freq':
            fig["layout"]['yaxis'].update(title="Modified <br> frequency")
        elif meth_traces.types[0] == 'bedgraph':
            fig["layout"]["yaxis"].update(title="Value")
        else:
            sys.exit(f"ERROR: unexpectedly not splitting for input of type {sample_type}")
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
    fig["layout"]["xaxis"].update(tickformat='g',
                                  separatethousands=True,
                                  range=[window.begin, window.end])
    fig["layout"].update(barmode='overlay',
                         title="Nucleotide modifications",
                         hovermode='closest',
                         plot_bgcolor='rgba(0,0,0,0)')
    if num_methrows > 10:
        for i in fig['layout']['annotations']:
            i['font']['size'] = 10
    utils.create_browser_output(fig, outfile, window)
    if static:
        import plotly.io as pio
        pio.write_image(fig, static, engine="kaleido")


if __name__ == '__main__':
    main()
