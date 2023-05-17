import astropy.units as u
import pytest
from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation

from rpnav.pulsar import Pulsar
from rpnav.observe.antenna import Antenna
from rpnav.observe.visualise import (
    plot_antenna_coverage,
    plot_flux_density,
    plot_pulsar_position,
)

@pytest.fixture(scope="module")
def parkes():
    return Antenna(
        name="Parkes telescope",
        signal_to_noise=10,  #
        system_temp=25.0,  # K
        astronomy_gain=0.6,  # (K/Jy)
        bandwidth=340e6,  # (Hz)
        centre_frequency=1374,  # (MHz)
        location=EarthLocation.from_geodetic(lat=-32.99327814611731, lon=148.26503125664433),
    )

@pytest.fixture(scope="module")
def fast():
    return Antenna(
        name="FAST Telescope",
        signal_to_noise=9.0,  #
        system_temp=25.0,  # (K)
        astronomy_gain=16.0,  # (K/Jy)
        bandwidth=512e6,  # (Hz)
        centre_frequency=1350,  # (MHz)
        location=EarthLocation.from_geodetic(lat=25.654006939684034, lon=106.85784898726897),
    )


@pytest.fixture(scope="module")
def kens():
    return Antenna(
        name="msfd_actual",
        signal_to_noise=1,
        system_temp=30.0,  # K
        bandwidth=40e6,  # (Hz)
        centre_frequency=1400,  # (MHz)
        location=EarthLocation.from_geodetic(lat=-33.77290643916046, lon=151.0976937264337),
        diameter=2,
    )

@pytest.fixture(scope="module")
def antenna(kens, parkes, fast):
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
                x=0.5
            )
        )
        plot_antenna_coverage(fig, antenna, t)
        # fig_pos.show()
        fig.show()


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
                x=0.5
            )
        )
        fig.show()


def test_return_az_el_visible(antenna):
    pulsars = Pulsar.load_catalogue(antenna.centre_frequency)
    t = Time.now()
    integ_time = 15 * 60
    pulsars.sort(key=lambda p: p.flux_density, reverse=True)
    for pulsar in pulsars:
        if pulsar.is_observable(antenna, t, integ_time):
            altaz = pulsar.get_alt_az(antenna, t)
            print(pulsar.name, pulsar.flux_density, altaz.alt, altaz.az)
