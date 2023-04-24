import astropy.units as u
import pandas as pd
import plotly.graph_objects as go
from astroplan import Observer
from astropy.coordinates import ICRS, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time, TimeDelta

from antenna import Antenna
from psrcat import calculate_visible_pulsars, load_catalogue


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
    zen = SkyCoord(
        90 * u.deg, 0 * u.deg, frame=AltAz(location=antenna.location, obstime=t)
    ).transform_to(ICRS)
    ra = float(zen.ra / (1.0 * u.deg))
    dec = float(zen.dec / (1.0 * u.deg))
    n_horiz = 90 - dec if dec > 0 else -90 + dec
    s_horiz = -90 + dec if dec > 0 else 90 - dec
    e_horiz = ra - 90 if ra > 90 else 360 - ra
    w_horiz = -90 + ra if ra > 270 else ra + 90
    fig.add_hrect(
        y0=s_horiz,
        y1=n_horiz,
        fillcolor="LightSalmon",
        opacity=0.3,
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


if __name__ == "__main__":
    """ """
    from config import fast, kens, parkes

    antenna = kens

    cat, flux_cat = load_catalogue(antenna.centre_freq)
    t = Time.now()
    t1 = t + TimeDelta(6.0 * 60 * 60 * u.s)
    for integ_time in [
        60,
    ]:  # 5*60, 15*60]:
        fig_dens = plot_flux_density(cat, flux_cat)
        vis = calculate_visible_pulsars(cat, flux_cat, antenna, t, integ_time)
        col = ["rgba(238, 118,0, 0.7)" if psr else "rgba(104,131,139, 0.4)" for psr in vis]
        fig_dens.add_hline(
            y=kens.min_observable_flux_density(integ_time),
            line_dash="dash",
            label=dict(text=antenna.name, textposition="end"),
        )
        fig_dens.update_traces(marker=dict(color=col))
        fig_dens.show()

        fig_pos = plot_pulsar_position(cat)
        fig_pos.update_traces(marker=dict(color=col))
        plot_antenna_coverage(fig_pos, antenna, t)
        plot_antenna_coverage(fig_pos, antenna, t1)
        fig_pos.show()
