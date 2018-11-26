import pandas as pd
import plotly
import plotly.graph_objs as go


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


def get_data():
    meth = pd.read_csv("methylation.txt", sep="\t")
    exp = [line.strip() for line in open("expanded-read-ids.txt")]
    nonexp = [line.strip() for line in open("non-expanded-read-ids.txt")]

    meth["allele"] = ""
    meth.loc[meth.read_name.isin(exp), "allele"] = "expanded"
    meth.loc[meth.read_name.isin(nonexp), "allele"] = "reference"
    return meth


meth = get_data()

html = plotly.offline.plot(
    {"data": [plot_read(meth, read) for read in meth.read_name.unique()],
     "layout": go.Layout(barmode='overlay',
                         title="Methylation per allele across the C9ORF72 locus")
     },
    output_type="div",
    show_link=False)

with open("plot.html", 'w') as output:
    output.write(html)
