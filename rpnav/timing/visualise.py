import plotly.graph_objects as go


def plot_residuals(resids: list[float], mjd: list[float], t_err: list[float]):
    """
    Plot timing residuals 

    Args:
        resids:     List of timing residuals
        mjd:        List of timestamps in format mean julian day
        t_err:      List of timing error

    Returns:
        Plotted residuals
    """
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


def plot_residuals_mse(x: list[float], y: list[float], residuals_rmse: list[float]):
    fig = go.Figure(data=[go.Surface(x=x, y=y, z=residuals_rmse)])
    fig.update_traces(
        contours_z=dict(show=True, usecolormap=True, highlightcolor="limegreen", project_z=True)
    )
    fig.update_layout(
        title="Timing Residuals RMSE Surface",
        autosize=False,
        scene_camera_eye=dict(x=1.87, y=0.88, z=-0.64),
        width=500,
        height=500,
        margin=dict(l=65, r=50, b=65, t=90),
    )
    return fig
