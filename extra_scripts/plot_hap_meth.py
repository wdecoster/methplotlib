from argparse import ArgumentParser
import pandas as pd
import plotly.graph_objects as go
from sklearn.decomposition import PCA
import re


def main():
    args = get_args()
    df = pd.read_csv(args.meth, sep="\t", index_col=0)
    df.columns = [re.sub('_v4.*', '', d) for d in df.columns]
    sample_info = pd.read_csv(args.sampleinfo, sep="\t", index_col=0)
    bars = bar_chart(df)
    pca = plot_pca(df, sample_info)
    with open("haplotype_specific_meth_plots.html", 'w') as output:
        output.write(bars)
        output.write(pca)


def bar_chart(df):
    uniq = df[df.count(axis=1) == 1]
    shared = df[df.count(axis=1) != 1]
    return go.Figure([go.Bar(x=shared.columns, y=shared.count(), name="Shared loci"),
                      go.Bar(x=uniq.columns, y=uniq.count(), name="Unique loci")],
                     layout=dict(title_text='Number of haplotype-specific loci', barmode='stack')) \
        .to_html(full_html=False, include_plotlyjs='cdn')


def plot_pca(df, sample_info):
    colors = color_per_sample(df.columns, sample_info)
    sklearn_pca = PCA(n_components=2)
    df = (df > 1).astype(int)
    pca = sklearn_pca.fit_transform(df.fillna(0).transpose())
    return go.Figure(data=[dict(type='scatter',
                                x=[pca[index, 0]],
                                y=[pca[index, 1]],
                                mode='markers',
                                name=name,
                                hoverinfo='name',
                                marker_color=color,
                                marker=dict(size=12, opacity=0.8))
                           for index, (name, color) in enumerate(zip(df.columns, colors))],
                     layout=dict(xaxis=dict(title='PC1', showline=False),
                                 yaxis=dict(title='PC2', showline=False),
                                 title='Principal Component Analysis after converting to binary')) \
        .to_html(full_html=False, include_plotlyjs='cdn')


def color_per_sample(columns, sample_info):
    from plotly.colors import DEFAULT_PLOTLY_COLORS as plcolors
    colors = dict(zip(sample_info['group'].unique(), plcolors))
    return [colors[sample_info.loc[d, 'group']] for d in columns]


def get_args():
    parser = ArgumentParser()
    parser.add_argument("meth", help="merged haplotype specific methylation file")
    parser.add_argument("sampleinfo", help="group information for included samples")
    return parser.parse_args()


if __name__ == '__main__':
    main()
