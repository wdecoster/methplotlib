from sklearn.decomposition import PCA
import plotly.express as px
import pandas as pd
import plotly
import plotly.graph_objs as go


def qc_plots(meth_data, window, qcpath=None, outpath=None):

    if qcpath is None and outpath is None:
        outfile = f"qc_report_{window.string}.html"
    elif qcpath is None:
        from pathlib import Path, PosixPath
        p = Path(outpath.format(region=window.string))
        Path.mkdir(p.parent, exist_ok=True, parents=True)
        outfile = str(p.parent / PosixPath("qc_" + p.stem + ".html"))
    else:
        from pathlib import Path
        p = Path(qcpath)
        Path.mkdir(p.parent, exist_ok=True, parents=True)
        outfile = qcpath

    with open(outfile, 'w') as qc_report:
        qc_report.write(num_sites_bar(meth_data))
        if len([m for m in meth_data if m.data_type == "nanopolish_freq"]) > 0:
            data = [m.table.rename({"methylated_frequency": m.name}, axis='columns')
                    for m in meth_data if m.data_type == "nanopolish_freq"]
            full = data[0].join(data[1:]).dropna(how="any", axis="index")
            qc_report.write(modified_fraction_histogram(full))
        if len([m for m in meth_data if m.data_type == "nanopolish_freq"]) > 2:
            qc_report.write(pairwise_correlation_plot(full))
            qc_report.write(pca(full))
            qc_report.write(global_box(data))
        if len([m for m in meth_data
                if m.data_type in ["nanopolish_call", "nanopolish_phased"]]) > 2:
            pass


def num_sites_bar(meth_data):
    trace = go.Bar(x=[m.name for m in meth_data],
                   y=[m.called_sites for m in meth_data])
    layout = dict(title="Number of called positions")
    return plotly.offline.plot(dict(data=[trace], layout=layout),
                               output_type="div",
                               show_link=False,
                               include_plotlyjs='cdn')


def pairwise_correlation_plot(full):
    trace = go.Splom(dimensions=[dict(label=l, values=full[l]) for l in full.columns],
                     marker=dict(size=4,
                                 line=dict(width=0.5,
                                           color='rgb(230,230,230)')),
                     diagonal=dict(visible=False)
                     )

    layout = go.Layout(
        title='Correlation of modification frequency',
        dragmode='select',
        width=1200,
        height=1200,
        autosize=False,
        hovermode=False,
        plot_bgcolor='rgba(240,240,240, 0.95)')

    for i in range(1, len(full.columns) + 1):
        layout[f"xaxis{i}"] = dict(showline=True, zeroline=False, gridcolor='#fff', ticklen=4)
        layout[f"yaxis{i}"] = dict(showline=True, zeroline=False, gridcolor='#fff', ticklen=4)

    return plotly.offline.plot(dict(data=[trace], layout=layout),
                               output_type="div",
                               show_link=False,
                               include_plotlyjs='cdn')


def pca(full):
    sklearn_pca = PCA(n_components=2)
    pca = sklearn_pca.fit_transform(full.transpose())
    data = [dict(type='scatter',
                 x=[pca[index, 0]],
                 y=[pca[index, 1]],
                 mode='markers',
                 name=name,
                 hoverinfo='name',
                 marker=dict(
                     size=12,
                     line=dict(
                         color='rgba(217, 217, 217, 0.14)',
                         width=0.5),
                     opacity=0.8))
            for index, name in enumerate(full.columns)]

    layout = dict(xaxis=dict(title='PC1', showline=False),
                  yaxis=dict(title='PC2', showline=False),
                  title="Principal Component Analysis"
                  )
    return plotly.offline.plot(dict(data=data, layout=layout),
                               output_type="div",
                               show_link=False,
                               include_plotlyjs='cdn')


def global_box(data):
    fig = px.box(pd.concat([d.reset_index(drop=True)
                            .rename({d.columns[0]: "freq"}, axis="columns")
                            .assign(dataset=d.columns[0]) for d in data], ignore_index=True),
                 x="dataset", y="freq", title="Global frequency of modification")
    return plotly.offline.plot(fig,
                               output_type="div",
                               show_link=False,
                               include_plotlyjs='cdn')


def modified_fraction_histogram(full):
    traces = [go.Histogram(x=full[dataset].dropna(),
                           histnorm='probability density',
                           xbins=dict(start=0, size=0.01, end=1),
                           name=dataset,
                           opacity=0.6
                           )
              for dataset in full.columns]
    layout = dict(barmode="overlay",
                  title="Histogram of modified fractions",
                  xaxis=dict(title="Modified fraction"),
                  yaxis=dict(title="Frequency"))
    return plotly.offline.plot(dict(data=traces,
                                    layout=layout),
                               output_type="div",
                               show_link=False,
                               include_plotlyjs='cdn')
