import pandas as pd
import plotly.graph_objects as go
from psrcat import load_catalogue
from astroplan import Observer
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

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
    """
    https://learn.astropy.org/tutorials/2-Coordinates-Transforms
    http://astro.wsu.edu/worthey/astro/html/lec-celestial-sph.html
    https://www.rpi.edu/dept/phys/observatory/obsastro6.pdf
    """
    lat = float(antenna.location.lat  / (1.*u.deg))
    lon = float(antenna.location.lon  / (1.*u.deg))
    n_horiz =  90 - lat  if lat > 0 else -90 + lat
    s_horiz =  -90 + lat if lat > 0 else 90 - lat
    fig.add_hrect(
        y0=s_horiz,
        y1=n_horiz,
        fillcolor="LightSalmon",
        opacity=0.5,
        layer="below",
        line_width=0,
    )
    fig.add_shape(
        type="circle",
        x0=lon - 90,
        x1=lon + 90,
        y0=s_horiz,
        y1=n_horiz,
        opacity=0.3,
        fillcolor="orange",
        line_color="orange"
    )
    fig.add_trace(
        go.Scatter(
            x=[lat,],
            y=[lon,],
            mode="markers",
            text=[antenna.name,],
            marker=dict(color="orange"),
        )
    )


def calculate_visible_pulsars(cat: Table, antenna: Antenna, t: Time, integ_time: float) -> list[bool]:
    vis: list[bool] = []
    for i in range(len(cat)):
        sky_pos = SkyCoord(frame="galactic", l=cat["GL"][i], b=cat["GB"][i], unit="deg")
        is_vis = (
            antenna.target_is_up(t, sky_pos) and cat["S1400"][i] > parkes.min_observable_flux_density(integ_time)
        )
        vis.append(is_vis)
    return vis


if __name__ == "__main__":
    """
    """
    from config import parkes, fast
    cat = load_catalogue()

    t = Time.now()
    for integ_time in [60,]: # 5*60, 15*60]:
        fig_dens = plot_flux_density(cat)
        vis = calculate_visible_pulsars(cat, parkes, t, integ_time)
        opacity = [1 if psr else 0.2 for psr in vis]
        fig_dens.add_hline(
            y=parkes.min_observable_flux_density(integ_time),
            line_dash="dash",
            label=dict(
                text=parkes.name,
                textposition="end"
            )
        )
        fig_dens.update_traces(
            marker=dict(opacity=opacity)
        )
        fig_dens.show()

        fig_pos = plot_pulsar_position(cat)
        fig_pos.update_traces(
            marker=dict(opacity=opacity)
        )
        plot_antenna_coverage(fig_pos, parkes)
        fig_pos.show()
