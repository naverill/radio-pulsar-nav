import astropy.units as u
import pytest
from astropy.time import Time, TimeDelta

from rpnav.observe.antenna import Antenna
from rpnav.observe.visualise import (
    plot_antenna,
    plot_antenna_coverage,
    plot_flux_density,
    plot_pulsar_position,
)
from rpnav.observe.pulsar import Pulsar
from rpnav.observe.observatories import PARKES, FAST, MSFD, WOODCHESTER_WEAK, WOODCHESTER_STRONG

@pytest.fixture(scope="module")
def antenna():
    return WOODCHESTER_WEAK

def test_plot_pulsar_position(antenna):
    pulsars = Pulsar.load_catalogue()
    # pulsars = Pulsar.load_catalogue(antenna.centre_frequency.to(u.MHz))
    t = Time.now()
    for integ_time in [60 * 60 * u.s, 6 * 60 * 60 * u.s, 24 * 60 * 60 * u.s]:
        fig = plot_pulsar_position(antenna, pulsars)
        vis = [pulsar.is_observable(antenna, t, integ_time) for pulsar in pulsars]
        col = ["rgba(238,118,0,0.7)" if psr else "rgba(104,131,139,0.4)" for psr in vis]
        fig.update_traces(marker=dict(color=col))
        fig.update_layout(
            title=dict(
                text=f"Pulsar Positioning over {integ_time} integration time",
                xanchor="center",
                yanchor="top",
                x=0.5,
            )
        )
        plot_antenna_coverage(fig, antenna, t)
        plot_antenna(fig, antenna, t)
        fig.show()
        # fig.write_image(f"outputs/{antenna.name}_pos_{integ_time}_f{antenna.centre_frequency.to_value(u.MHz)}.png")


def test_plot_pulsar_horizon_integ_time(antenna):
    pulsars = Pulsar.load_catalogue(antenna.centre_frequency.to(u.MHz))
    t = Time.now()
    below_horiz = [pulsar.is_observable(antenna, t, 1e10 * u.s) for pulsar in pulsars]
    for integ_time in [5 * 60 * u.s, 15 * 60 * u.s, 60 * 60 * u.s]:
        fig = plot_pulsar_position(antenna, pulsars)
        vis = [pulsar.is_observable(antenna, t, integ_time) for pulsar in pulsars]
        col = []
        for is_obs, is_vis in zip(below_horiz, vis):
            if is_obs and is_vis:
                col.append("rgba(102,205,170,0.8)")
            elif is_obs:
                col.append("rgba(238,118,0,0.8)")
            else:
                col.append("rgba(104,131,139,0.8)")

        fig.update_traces(marker=dict(color=col))
        fig.update_layout(
            title=dict(
                text=f"Pulsar Positioning over {integ_time} integration time",
                xanchor="center",
                yanchor="top",
                x=0.5,
            )
        )
        plot_antenna(fig, antenna, t)
        fig.update_layout({"paper_bgcolor": "rgba(0,0,0,0)"})
        fig.update_xaxes(showgrid=True, gridwidth=1, minor_ticks="inside")
        fig.update_yaxes(showgrid=True, gridwidth=1, minor_ticks="inside")
        fig.show()
        # fig.write_image(f"outputs/{antenna.name}_horiz_{integ_time}_f{antenna.centre_frequency.to_value(u.MHz)}.png")


def test_plot_pulsar_horizon_obs_time(antenna):
    pulsars = Pulsar.load_catalogue(antenna.centre_frequency.to(u.MHz))
    integ_time = 15 * 60  * u.s
    for t_delta in [0  * u.s, 6 * 60 * 60  * u.s, 12 * 60 * 60  * u.s, 18 * 60 * 60  * u.s]:
        t = Time.now() - TimeDelta(t_delta)
        fig = plot_pulsar_position(antenna, pulsars)
        vis = [pulsar.is_observable(antenna, t, integ_time) for pulsar in pulsars]
        below_horiz = [pulsar.is_observable(antenna, t, 1e10  * u.s) for pulsar in pulsars]
        col = []
        for is_obs, is_vis in zip(below_horiz, vis):
            if is_obs and is_vis:
                col.append("rgba(102,205,170,0.8)")
            elif is_obs:
                col.append("rgba(238,118,0,0.8)")
            else:
                col.append("rgba(104,131,139,0.8)")

        fig.update_traces(marker=dict(color=col))
        fig.update_layout(
            title=dict(
                text=f"Pulsar Positioning over {integ_time} integration time",
                xanchor="center",
                yanchor="top",
                x=0.5,
            )
        )
        # plot_antenna(fig, antenna, t)
        fig.update_layout({"paper_bgcolor": "rgba(0,0,0,0)"})
        fig.update_xaxes(showgrid=True, gridwidth=1, minor_ticks="inside")
        fig.update_yaxes(showgrid=True, gridwidth=1, minor_ticks="inside")
        fig.show()
        # fig.write_image(f"outputs/{antenna.name}_obs_{t_delta}_f{antenna.centre_frequency.to_value(u.MHz)}.png")



def test_plot_pulsar_observability(antenna):
    pulsars = Pulsar.load_catalogue(antenna.centre_frequency.to(u.MHz))
    t = Time.now()

    for integ_time in [60 * 60 * u.s, 6 * 60 * 60 * u.s, 24 * 60 * 60 * u.s]:
        fig = plot_flux_density(antenna, pulsars)
        vis = [pulsar.is_observable(antenna, t, integ_time) for pulsar in pulsars]
        col = ["rgba(238,118,0,0.7)" if psr else "rgba(104,131,139,0.4)" for psr in vis]

        fig.add_hline(
            y=antenna.min_observable_flux_density(integ_time).to_value(u.mJy),
            line_dash="dash",
            label=dict(text=antenna.name, textposition="end"),
        )
        fig.update_traces(marker=dict(color=col))
        fig.update_layout(
            title=dict(
                text=f"Pulsar observability over {integ_time} integration time",
                xanchor="center",
                yanchor="top",
                x=0.5,
            )
        )
        # fig.show()
        # fig.write_image(f"outputs/{antenna.name}_int_{integ_time}_f{antenna.centre_frequency.to_value(u.MHz)}.png")


def test_return_az_el_visible(antenna):
    pulsars = Pulsar.load_catalogue(antenna.centre_frequency.to(u.MHz))
    t = Time.now()
    integ_time = 15 * 60 * u.s
    pulsars.sort(key=lambda p: p.flux_density, reverse=True)
    for pulsar in pulsars:
        if pulsar.is_observable(antenna, t, integ_time):
            altaz = pulsar.get_alt_az(antenna, t)
            print(pulsar.name, pulsar.flux_density, altaz.alt, altaz.az)

def test_get_toa(antenna):
    print(PARKES.snr())
    print(WOODCHESTER_WEAK.toa_err(PARKES, ref_integration_time = 1800 * u.s, integration_time=60 * 60 * u.s))
    print(WOODCHESTER_STRONG.toa_err(PARKES, ref_integration_time = 1800 * u.s, integration_time=60 * 60 * u.s))