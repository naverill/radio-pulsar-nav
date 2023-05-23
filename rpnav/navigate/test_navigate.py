from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from astropy.time import Time
from pint import logging
from pint.fitter import DownhillWLSFitter
from pint.models import get_model_and_toas

from rpnav.constants import SPEED_OF_LIGHT
from rpnav.navigate.gradient_descent import localise
from rpnav.observe.observer import Observer

logging.setup(level="ERROR")
FILE_DIR = Path(__file__).parent


class Args:
    def __init__(
        self, timfile: str, parfile: str, est: tuple[float], frame: str, err: tuple[float]
    ):
        self.timfile = timfile
        self.parfile = parfile
        self.est = est
        self.frame = frame
        self.err = err


def process_estimate(args) -> EarthLocation:
    if len(args.est) < 2 or len(args.est) > 3:
        logging.error("Position must contain estimates for x and y axes")

    if args.frame == "geodetic":
        if len(args.est == 2):
            est_loc = EarthLocation.from_geodetic(lat=args.est[0], lon=args.est[1])
        else:
            est_loc = EarthLocation.from_geodetic(
                lat=args.est[0], lon=args.est[1], height=args.est[2]
            )
    elif args.frame == "geocentric":
        est_loc = EarthLocation.from_geocentric(
            x=args.est[0] * u.m, y=args.est[1] * u.m, z=args.est[2] * u.m
        )
    else:
        logging.error("Invalid reference frame. Must be geodetic or geocentric")
    return est_loc


@pytest.fixture(scope="module")
def parkes():
    sim = ("parkes_sim",)
    return Args(
        timfile=f"{FILE_DIR}/../simulate/{sim}/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/../simulate/{sim}/output/real_0/J0835-4510.tdb.par",
        est=[-4554231.5, 2816759.1, -3454036.3],
        frame="geocentric",
        err=None,
    )


@pytest.fixture(scope="module")
def msfd():
    sim = "msfd_sim"
    loc = EarthLocation.from_geodetic(
        lat=-33.77290643916046 * u.deg, lon=151.0976937264337 * u.deg
    ).geocentric
    return Args(
        timfile=f"{FILE_DIR}/../simulate/{sim}/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/../simulate/{sim}/output/real_0/J0835-4510.tdb.par",
        est=[loc[0].to_value(u.m), loc[1].to_value(u.m), loc[2].to_value(u.m)],
        frame="geocentric",
        err=None,
    )


def test_gradient_descent(msfd):
    half_pulsar_width = (0.08932838502359318 * u.s / 2) * SPEED_OF_LIGHT
    true_loc = process_estimate(msfd)
    loc_err = EarthLocation.from_geocentric(
        x=true_loc.x + np.random.normal(scale=np.sqrt(half_pulsar_width.to_value(u.m) * 2)) * u.m,
        y=true_loc.y + np.random.normal(scale=np.sqrt(half_pulsar_width.to_value(u.m) * 2)) * u.m,
        z=true_loc.z,
    )
    start = (loc_err.geodetic.lon.to_value(u.deg), loc_err.geodetic.lat.to_value(u.deg))
    msfd_sim = Observer(
        name="msfd_actual",
        itoa_code="MSFD",
        location=true_loc,
        origin="True location of CSIRO Marsfield telescope site",
        time=Time(60002.3, format="mjd"),
    )
    msfd_sim.update()
    msfd_est = Observer(
        name="msfd_estimated",
        itoa_code="MSFD_EST",
        origin="Estimated location of CSIRO Marsfield site",
        location=loc_err,
        time=Time(60002.3, format="mjd"),
    )
    msfd_est.update()

    model, toas = get_model_and_toas(msfd.parfile, msfd.timfile)
    fitter = DownhillWLSFitter(toas, model)
    est_loc = localise(start, msfd_est, fitter)
    est_loc = EarthLocation.from_geodetic(lon=est_loc[0], lat=est_loc[1])
    print(
        np.abs(true_loc.x - est_loc.x),
        np.abs(true_loc.y - est_loc.y),
    )
