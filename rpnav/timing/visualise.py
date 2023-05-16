import plotly.graph_objects as go


def plot_residuals(resids: list[float], mjd: list[float], t_err: list[float]):
    print(len(resids), len(mjd), len(t_err))
    fig = go.Figure(
        data=go.Scatter(
            x=mjd,
            y=resids,
            error_y=dict(
                type="data",  # value of error bar given in data coordinates
                array=t_err,
                visible=True,
            ),
            mode="markers",
        )
    )
    return fig
