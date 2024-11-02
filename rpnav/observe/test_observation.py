import astropy.units as u
import pytest
from astropy.time import Time, TimeDelta

from rpnav.observe.observation import Observation
from rpnav.observe.pulsar import Pulsar
from rpnav.observe.observatories import PARKES, FAST, MSFD, WOODCHESTER
from rpnav.observe.antenna import Antenna
from astropy import constants as const


@pytest.fixture(scope="module")
def antenna():
    return WOODCHESTER

def test_get_toa(antenna):
    observations = []

    """
    Reference antenna defined in Feasibility Study for Pulsar Tim
    """
    SALA = Antenna(
        name="SALA",
        centre_frequency=1 * u.GHz,
        bandwidth=200 * u.MHz,
        effective_area=10 * u.m * u.m,
        location=PARKES.location,
        time=PARKES.time,
    )

    pulsars = Pulsar.load_catalogue()
    for p in pulsars:
        if p.name == "J1909-3744":
            observations.append(
                Observation(
                    PARKES,
                    p,
                    snr=146.09 * u.dimensionless_unscaled,
                    toa_err=(145 * u.ns).to(u.s),
                    integration_time= 1800 * u.s
                )
            )
        elif p.name == "B1937+21": 
            observations.append(
                Observation(
                    SALA,
                    p,
                    snr= (-55.6 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(3162 * u.m / const.c).to(u.s),
                    integration_time=(36.36 * u.min).to(u.s)
                )
            )
        elif p.name == "B0736-40": 
            observations.append(
                Observation(
                    SALA,
                    p,
                    snr= (-50.2 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(9246467 * u.m / const.c).to(u.s),
                    integration_time=(3.07 * u.min).to(u.s)
                )
            )
        elif p.name == "B1451-68":
            observations.append(
                Observation(
                    SALA,
                    p,
                    snr= (-50.0 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(452529 * u.m / const.c).to(u.s),
                    integration_time=(2.35 * u.min).to(u.s)
                )
            )
        elif p.name == "B0950-08": 
            observations.append(
                Observation(
                    SALA,
                    p,
                    snr= (-49.8 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(316174 * u.m / const.c).to(u.s),
                    integration_time=(1.04 * u.min).to(u.s)
                )
            )
        elif p.name == "B0329+54": 
            observations.append(
                Observation(
                    observer=SALA,
                    pulsar=p,
                    snr= (-45.2 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(290556 * u.m / const.c).to(u.s),
                    integration_time=(0.12 * u.min).to(u.s)
                )
            )

    for obs in observations:
        wood_weak = Observation(
            observer=WOODCHESTER,
            pulsar=obs.pulsar,
            snr=5 * u.dimensionless_unscaled,
            integration_time=60 * 60 * u.s
        )
        wood_strong = Observation(
            observer=WOODCHESTER,
            pulsar=obs.pulsar,
            snr=20 * u.dimensionless_unscaled,
            integration_time=60 * 60 * u.s
        )

        print(wood_weak.toa_err(reference=obs))
        print(wood_strong.toa_err(reference=obs))