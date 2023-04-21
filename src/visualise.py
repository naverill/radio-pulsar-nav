import pandas as pd
import plotly.graph_objects as go
from psrcat import load_catalogue
from scipy.special import jv

from antenna import Antenna

def plot_flux_density(pulsar_cat: pd.DataFrame) -> go.Figure:
    fig = go.Figure(
        data = go.Scatter(
            x=pulsar_cat["P0"],
            y=pulsar_cat["S1400"],
            mode="markers",
            text=pulsar_cat["NAME"]
        )
    )
    fig.update_layout(
        title="Pulsar observability",
        xaxis_title="Pular period (s)",
        yaxis_title="Mean Flux Density (mJ)",
        xaxis_type = "log",
        yaxis_type = "log"
    )
    return fig

def plot_pulsar_position(pulsar_cat: pd.DataFrame) -> go.Figure:
    fig = go.Figure(
        data = go.Scatter(
            x=pulsar_cat["RAJD"],
            y=pulsar_cat["DECJD"],
            mode="markers",
            text=pulsar_cat["NAME"]
        )
    )
    fig.update_layout(
        title="Pulsar Positioning",
        xaxis_title="Right ascension (J2000) (degrees)",
        yaxis_title="Declination (J2000) (degrees)",
        yaxis_range=[-180, 180],
        xaxis_range=[0, 360],
    )
    return fig

def plot_antenna_coverage(fig: go.Figure, antenna: Antenna):
    fig.add_trace(
        go.Scatter(
            x=[antenna.lat,],
            y=[antenna.long,],
            mode="markers",
            text=[antenna.name,],
            marker=dict(color="orange"),
        )
    )
    fig.add_hrect(
        y0=antenna.long - antenna.max_abs_galactic_lat,
        y1=antenna.long + antenna.max_abs_galactic_lat,
        fillcolor="LightSalmon",
        opacity=0.5,
        layer="below",
        line_width=0,
    )
    fig.add_shape(
        type="circle",
        x0=antenna.lat - antenna.max_abs_galactic_lat,
        x1=antenna.lat + antenna.max_abs_galactic_lat,
        y0=antenna.long - antenna.max_abs_galactic_lat,
        y1=antenna.long + antenna.max_abs_galactic_lat,
        opacity=0.3,
        fillcolor="orange",
        line_color="orange"
    )

if __name__ == "__main__":
    from config import parkes, fast
    cat = load_catalogue(parkes)
    fig = plot_pulsar_position(cat)
    plot_antenna_coverage(fig, parkes)
    fig.show()

    fig = plot_flux_density(cat)
    for integ_time in [15*60,]: #10, 60, 600, 3600]:
        fig.add_hline(
            y=parkes.min_observable_flux_density(integ_time),
            line_dash="dash",
            label=dict(
                text=parkes.name,
                textposition="end"
            )
        )
        fig.add_hline(
            y=fast.min_observable_flux_density(integ_time),
            line_dash="dash",
            label=dict(
                text=fast.name,
                textposition="end"
            )
        )
    fig.show()
