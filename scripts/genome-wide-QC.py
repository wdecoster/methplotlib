from methplotlib.import_methylation import get_data
from argparse import ArgumentParser
import sys
import methplotlib.qc as qc


def main():
    args = get_args()
    meth_data = get_data(args.methylation, args.names, None, 5)
    data = [m.table.rename({"methylated_frequency": m.name}, axis='columns')
            for m in meth_data if m.data_type == "frequency"]
    labels = [m.name for m in meth_data if m.data_type == "frequency"]
    full = data[0].join(data[1:]).dropna(how="any", axis="index")
    with open("qc_report_genome-wide.html", 'w') as qc_report:
        qc_report.write(qc.num_sites_bar(meth_data))
        qc_report.write(qc.pairwise_correlation_plot(full, labels))
        qc_report.write(qc.pca(full, labels))
        qc_report.write(qc.global_box(data))


def get_args():
    parser = ArgumentParser(description="genome wide splom")
    parser.add_argument("-m", "--methylation",
                        nargs='+',
                        help="nanopolish methylation calls or frequency output",
                        required=True)
    parser.add_argument("-n", "--names",
                        nargs='+',
                        help="names of datasets in --methylation",
                        required=True)
    args = parser.parse_args()
    if not len(args.names) == len(args.methylation):
        sys.exit("INPUT ERROR: Expecting the same number of names as datasets!")
    return args


if __name__ == '__main__':
    main()
