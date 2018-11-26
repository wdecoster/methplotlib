import pandas as pd
import plotly
import plotly.graph_objs as go
from argparse import ArgumentParser


def plot_read(df, read):
    if df.loc[df.read_name == read, "allele"].unique() == "expanded":
        col = ('rgb(255,0,0)')
    else:
        col = ('rgb(0,0,255)')
    return go.Scatter(
        x=df.loc[df.read_name == read, "start"],
        y=df.loc[df.read_name == read, "log_lik_ratio"],
        mode='lines+markers',
        name=read,
        line=dict(color=col)
    )


def get_data(methylation, list_file1, list_file2):
    meth = pd.read_csv(methylation, sep="\t")
    exp = [line.strip() for line in open(list_file1)]
    nonexp = [line.strip() for line in open(list_file2)]

    meth["allele"] = ""
    meth.loc[meth.read_name.isin(exp), "allele"] = "expanded"
    meth.loc[meth.read_name.isin(nonexp), "allele"] = "reference"
    return meth


def get_args():
    parser = ArgumentParser(
        description="plot per read likelihood of methylation split by two files")
    parser.add_argument("methylation", help="output of nanopolish call methylation")
    parser.add_argument("--list1", help="list of readIDs for group 1", required=True)
    parser.add_argument("--list2", help="list of readIDs for group 2", required=True)
    parser.add_argument("--output", help="name of output html file",
                        default="per_read_methylation_plot.html")
    parser.add_argument("--legend", help="show legend", action="store_true")
    return parser.parse_args()


def main():
    args = get_args()
    meth = get_data(args.methylation, args.list1, args.list2)

    html = plotly.offline.plot(
        {"data": [plot_read(meth, read) for read in meth.read_name.unique()],
         "layout": go.Layout(barmode='overlay',
                             title="Methylation per allele across the C9ORF72 locus",
                             xaxis=dict(title='Genomic position'),
                             yaxis=dict(title='Log likelihood of methylated CpG'),
                             showLegend=args.legend,
                             )
         },
        output_type="div",
        show_link=False)

    with open(args.output, 'w') as output:
        output.write(html)


if __name__ == '__main__':
    main()
