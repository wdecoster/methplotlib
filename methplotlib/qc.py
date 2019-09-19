from sklearn.decomposition import PCA
import plotly.express as px
import pandas as pd
import plotly
import plotly.graph_objs as go


def pairwise_correlation_plot(full, labels, window):
    trace = go.Splom(dimensions=[dict(label=l, values=full[l]) for l in labels],
                     marker=dict(size=4,
                                 line=dict(width=0.5,
                                           color='rgb(230,230,230)')),
                     diagonal=dict(visible=False)
                     )

    layout = go.Layout(
        title='Correlation of methylation frequency',
        dragmode='select',
        width=1200,
        height=1200,
        autosize=False,
        hovermode=False,
        plot_bgcolor='rgba(240,240,240, 0.95)')

    for i in range(1, len(full.columns) + 1):
        layout["xaxis{}".format(i)] = dict(
            showline=True, zeroline=False, gridcolor='#fff', ticklen=4)
        layout["yaxis{}".format(i)] = dict(
            showline=True, zeroline=False, gridcolor='#fff', ticklen=4)

    with open("methylation_frequency_correlation_{}.html".format(window.string), 'w') as output:
        output.write(plotly.offline.plot(dict(data=[trace], layout=layout),
                                         output_type="div",
                                         show_link=False)
                     )


def pca(full, labels, window):
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
            for index, name in enumerate(labels)]

    layout = dict(
        xaxis=dict(title='PC1', showline=False),
        yaxis=dict(title='PC2', showline=False)
    )
    with open("principal_components_{}.html".format(window.string), 'w') as output:
        output.write(plotly.offline.plot(dict(data=data, layout=layout),
                                         output_type="div",
                                         show_link=False)
                     )


def global_box(data, window):
    fig = px.box(pd.concat([d.reset_index(drop=True)
                            .rename({d.columns[0]: "freq"}, axis="columns")
                            .assign(dataset=d.columns[0]) for d in data], ignore_index=True),
                 y="freq", x="dataset")
    with open("global_box_plot_{}.html".format(window.string), 'w') as output:
        output.write(plotly.offline.plot(fig,
                                         output_type="div",
                                         show_link=False)
                     )
