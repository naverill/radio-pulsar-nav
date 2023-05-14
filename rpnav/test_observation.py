import astropy.units as u
import pytest
from astropy.time import Time, TimeDelta

from rpnav.config import fast, kens, parkes
from rpnav.pulsar import Pulsar
from rpnav.visualise import (
    plot_antenna_coverage,
    plot_flux_density,
    plot_pulsar_position,
)


@pytest.fixture(scope="module")
def antenna():
    return kens


def test_plot_pulsar_position(antenna):
    pulsars = Pulsar.load_catalogue(antenna.centre_frequency)
    t = Time.now()
    for integ_time in [60 * 60, 6 * 60 * 60, 24 * 60 * 60]:
        fig = plot_pulsar_position(pulsars)
        vis = [pulsar.is_observable(antenna, t, integ_time) for pulsar in pulsars]
        col = ["rgba(238,118,0,0.7)" if psr else "rgba(104,131,139,0.4)" for psr in vis]
        fig.update_traces(marker=dict(color=col))
        fig.update_layout(
            title=dict(
                text=f"Pulsar Positioning over {integ_time}s integration time",
                xanchor="center",
                yanchor="top",
            )
        )
        plot_antenna_coverage(fig, antenna, t)
        # fig_pos.show()
        fig.write_image("tmp.png")


def test_plot_pulsar_observability(antenna):
    pulsars = Pulsar.load_catalogue(antenna.centre_frequency)
    t = Time.now()
    for integ_time in [60 * 60, 6 * 60 * 60, 24 * 60 * 60]:
        fig = plot_flux_density(pulsars)
        vis = [pulsar.is_observable(antenna, t, integ_time) for pulsar in pulsars]
        col = ["rgba(238,118,0,0.7)" if psr else "rgba(104,131,139,0.4)" for psr in vis]
        fig.add_hline(
            y=antenna.min_observable_flux_density(integ_time),
            line_dash="dash",
            label=dict(text=antenna.name, textposition="end"),
        )
        fig.update_traces(marker=dict(color=col))
        fig.update_layout(
            title=dict(
                text=f"Pulsar observability over {integ_time}s integration time",
                xanchor="center",
                yanchor="top",
            )
        )
        # fig.show()
        fig.write_image("tmp.png")


def test_return_az_el_visible(antenna):
    pulsars = Pulsar.load_catalogue(antenna.centre_frequency)
    t = Time.now()
    integ_time = 15 * 60
    pulsars.sort(key=lambda p: p.flux_density, reverse=True)
    for pulsar in pulsars:
        if pulsar.is_observable(antenna, t, integ_time):
            altaz = pulsar.get_alt_az(antenna, t)
            print(pulsar.name, pulsar.flux_density, altaz.alt, altaz.az)
