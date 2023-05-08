import astropy.units as u
import pandas as pd
import plotly.graph_objects as go
from astroplan import Observer
from astropy.coordinates import ICRS, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time

from radio_pulsar_nav.antenna import Antenna


def plot_flux_density(pulsar_cat: pd.DataFrame, flux_cat) -> go.Figure:
    fig = go.Figure(
        data=go.Scatter(
            x=pulsar_cat["P0"], y=flux_cat["S"], mode="markers", text=pulsar_cat["NAME"]
        )
    )
    fig.update_layout(
        title="Pulsar observability",
        xaxis_title="Pular period (s)",
        yaxis_title="Mean Flux Density (mJ)",
        xaxis_type="log",
        yaxis_type="log",
    )
    return fig


def plot_pulsar_position(pulsar_cat: pd.DataFrame) -> go.Figure:
    fig = go.Figure(
        data=go.Scatter(
            x=pulsar_cat["RAJD"],
            y=pulsar_cat["DECJD"],
            mode="markers",
            text=pulsar_cat["NAME"],
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


def plot_antenna_coverage(fig: go.Figure, antenna: Antenna, t: Time):
    """
    https://learn.astropy.org/tutorials/2-Coordinates-Transforms
    http://astro.wsu.edu/worthey/astro/html/lec-celestial-sph.html
    https://www.rpi.edu/dept/phys/observatory/obsastro6.pdf
    """
    zen = antenna.location.get_itrs(t).transform_to(ICRS)
    ra = float(zen.ra / (1.0 * u.deg))
    dec = float(zen.dec / (1.0 * u.deg))
    n_horiz = (dec + 90 + 180) % 360 - 180
    s_horiz = (dec - 90 + 180) % 360 - 180
    w_horiz = (ra + 90) % 360
    e_horiz = (ra - 90) % 360
    fig.add_hrect(
        y0=s_horiz,
        y1=n_horiz,
        fillcolor="LightSalmon",
        opacity=0.2,
        layer="below",
        line_width=0,
    )
    fig.add_shape(
        type="circle",
        x0=e_horiz,
        x1=w_horiz,
        y0=s_horiz,
        y1=n_horiz,
        opacity=0.1,
        fillcolor="orange",
        line_color="orange",
    )
    fig.add_trace(
        go.Scatter(
            x=[
                ra,
            ],
            y=[
                dec,
            ],
            mode="markers",
            text=[
                antenna.name,
            ],
            marker=dict(color="midnightblue"),
        )
    )
