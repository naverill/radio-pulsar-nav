import astropy.units as u
import pytest
from astropy.time import Time, TimeDelta

from rpnav.observe.observation import Observation
from rpnav.observe.pulsar import Pulsar
from rpnav.observe.observatories import PARKES, FAST, MSFD, WOODCHESTER

@pytest.fixture(scope="module")
def antenna():
    return WOODCHESTER

def test_get_toa(antenna):
    psr = None
    pulsars = Pulsar.load_catalogue()
    for p in psr:
        if p.name == "J1909–3744":
            psr = p

    if psr is None:
        return
    
    reference = Observation(
        PARKES,
        psr,
        snr=10 * u.dimensionless_unscaled,
        toa_err=145 * u.ns,
        integration_time= 1800 * u.s
    )
    wood_weak = Observation(
        WOODCHESTER,
        psr,
        snr=5 * u.dimensionless_unscaled,
        integration_time=60 * 60 * u.s
    )
    wood_strong = Observation(
        WOODCHESTER,
        psr,
        snr=20 * u.dimensionless_unscaled,
        integration_time=60 * 60 * u.s
    )

    print(wood_weak.toa_err(reference=reference))
    print(wood_strong.toa_err(reference=reference))