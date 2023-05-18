from math import log2

import astropy.units as u
import pandas as pd
import plotly.graph_objects as go
from astropy.coordinates import ICRS, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time

from rpnav.observe.antenna import Antenna
from rpnav.pulsar import Pulsar


def plot_flux_density(pulsars: list[Pulsar]) -> go.Figure:
    widths: list[float] = []
    fluxs: list[float] = []
    names: list[float] = []
    for p in pulsars:
        widths.append(p.pulse_width_10)
        fluxs.append(p.flux_density)
        names.append(p.name)

    fig = go.Figure(data=go.Scatter(x=widths, y=fluxs, mode="markers", text=names))
    fig.update_layout(
        title="Pulsar observability",
        xaxis_title="Pular period (s)",
        yaxis_title="Mean Flux Density (mJ)",
        xaxis_type="log",
        yaxis_type="log",
    )
    return fig


def plot_pulsar_position(pulsars: list[Pulsar]) -> go.Figure:
    ras: list[float] = []
    decs: list[float] = []
    names: list[float] = []
    size: list[float] = []
    flux_dens: list[float] = []
    for p in pulsars:
        ras.append(p.ra.value)
        decs.append(p.dec.value)
        names.append(p.name)
        flux_dens.append(p.flux_density)
        size.append(log2(p.flux_density) * 2 if log2(p.flux_density) * 2 > 2 else 2)

    fig = go.Figure(data=go.Scatter(x=ras, y=decs, mode="markers", text=names, marker_size=size))
    fig.update_layout(
        title="Pulsar Positioning",
        xaxis_title="Right ascension (J2000) (degrees)",
        yaxis_title="Declination (J2000) (degrees)",
        yaxis_range=[-90, 90],
        xaxis_range=[0, 360],
    )
    return fig


def plot_antenna(fig: go.Figure, antenna: Antenna, t: Time):
    zen = antenna.location.get_itrs(t).transform_to(ICRS)
    ra = zen.ra.value
    dec = zen.dec.value
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
            marker=dict(color="rgba(69, 139, 116, 1)"),
        )
    )
    return fig


def plot_antenna_coverage(fig: go.Figure, antenna: Antenna, t: Time):
    """
    https://learn.astropy.org/tutorials/2-Coordinates-Transforms
    http://astro.wsu.edu/worthey/astro/html/lec-celestial-sph.html
    https://www.rpi.edu/dept/phys/observatory/obsastro6.pdf
    """
    zen = antenna.location.get_itrs(t).transform_to(ICRS)
    ra = zen.ra.value
    dec = zen.dec.value
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
    return fig
