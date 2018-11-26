import pandas as pd
import plotly
import plotly.graph_objs as go
from argparse import ArgumentParser


def plot_read(df, read):
    if df.loc[df.read_name == read, "allele"].unique() == "a1":
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


def get_data(args):
    meth = pd.read_csv(args.methylation, sep="\t")
    a1 = [line.strip() for line in open(args.allele1)]
    a2 = [line.strip() for line in open(args.allele2)]

    meth["allele"] = ""
    meth.loc[meth.read_name.isin(a1), "allele"] = "a1"
    meth.loc[meth.read_name.isin(a2), "allele"] = "a2"
    return meth.loc[meth["allele"] != ""]


def single_observation_plot(meth):
    html = plotly.offline.plot(
        {"data": [plot_read(meth, read) for read in meth.read_name.unique()],
         "layout": go.Layout(barmode='overlay',
                             title="Methylation per allele")
         },
        output_type="div",
        show_link=False)

    with open("single_observation_plot.html", 'w') as output:
        output.write(html)


def windowed_allele(meth, allele):
    return meth.loc[meth.allele == allele, ['log_lik_ratio', 'start']] \
        .sort_values('start') \
        .groupby('start') \
        .mean() \
        .rolling(window=10, center=True) \
        .mean()


def windowed_mean_plot(meth):
    alleles = [windowed_allele(meth, "a1"), windowed_allele(meth, "a2")]
    colors = [('rgb(0,0,255)'), ('rgb(255,0,0)')]
    data = [go.Scatter(x=a.index, y=a.loc[:, "log_lik_ratio"], mode='lines+markers',
                       line=dict(color=c)) for a, c in zip(alleles, colors)]

    html = plotly.offline.plot(
        {"data": data,
         "layout": go.Layout(barmode='overlay',
                             title="Methylation per allele")
         },
        output_type="div",
        show_link=False)

    with open("windowed_mean_plot.html", 'w') as output:
        output.write(html)


def main():
    meth = get_data(get_args())
    single_observation_plot(meth)
    windowed_mean_plot(meth)


def get_args():
    parser = ArgumentParser(description="plotting methylation for two alleles separately")
    parser.add_argument("methylation", help="output of nanopolish call-methylation")
    parser.add_argument("--allele1", required=True, help="readIds of allele1")
    parser.add_argument("--allele2", required=True, help="readIds of allele2")
    return parser.parse_args()


if __name__ == '__main__':
    main()
