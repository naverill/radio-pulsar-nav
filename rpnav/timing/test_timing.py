import argparse
import json
from math import pi
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from astropy.time import Time
from pint.observatory import get_observatory

from rpnav.constants import SPEED_OF_LIGHT
from rpnav.observe.antenna import Antenna
from rpnav.observe.observer import Observer
from rpnav.timing.timing import fit_residuals
from rpnav.timing.visualise import plot_residuals, plot_residuals_mse

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
        raise parser.error("Position must contain estimates for x and y axes")

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
        parser.error("Invalid reference frame. Must be geodetic or geocentric")
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


def test_residuals_parkes_correct_position(parkes):
    est_loc = process_estimate(parkes)
    pks_sim = Observer(
        name="pks",
        itoa_code="PARKES",
        location=est_loc,
        origin="Actual location of Parkes telescope site",
        time=Time(60002.3, format="mjd"),
    )
    pks_sim.update()
    fitter = fit_residuals(parkes.parfile, parkes.timfile, parkes.est, parkes.err)
    fig = plot_residuals(
        fitter.resids.time_resids.to_value(u.us).astype(float),
        fitter.toas.get_mjds().value.astype(float),
        fitter.toas.get_errors().to_value(u.us).astype(float),
    )
    fig.update_layout(width=1000, height=1000)
    fig.write_image("parkes.png")
    # fig.show()


def test_residuals_msfd_true(msfd):
    est_loc = process_estimate(msfd)
    msfd_sim = Observer(
        name="msfd_actual",
        itoa_code="MSFD",
        location=EarthLocation.from_geocentric(
            x=-4646251.90 * u.m, y=2565112.87 * u.m, z=-3525535.83 * u.m
        ),
        origin="Actual location of Parkes telescope site",
        time=Time(60002.3, format="mjd"),
    )
    msfd_sim.update()
    msfd_est = Observer(
        name="msfd_estimated",
        itoa_code="MSFD_EST",
        origin="Estimated location of CSIRO Marsfield site",
        location=est_loc,
        time=Time(60002.3, format="mjd"),
    )
    msfd_est.update()
    fitter = fit_residuals(msfd.parfile, msfd.timfile, msfd.est, msfd.err)
    fig = plot_residuals(
        fitter.resids.time_resids.to_value(u.us).astype(float),
        fitter.toas.get_mjds().value.astype(float),
        fitter.toas.get_errors().to_value(u.us).astype(float),
    )
    fig.update_layout(width=1000, height=1000)
    fig.write_image("msfd.png")
    # fig.show()


def test_residuals_msfd_incorrect(msfd):
    msfd.est = (msfd.est[0], msfd.est[1], msfd.est[2] + 10000)
    est_loc = process_estimate(msfd)
    msfd_sim = Observer(
        name="msfd_actual",
        itoa_code="MSFD",
        location=EarthLocation.from_geocentric(
            x=-4646251.90 * u.m, y=2565112.87 * u.m, z=-3525535.83 * u.m
        ),
        origin="True location of CSIRO Marsfield telescope site",
        time=Time(60002.3, format="mjd"),
    )
    msfd_sim.update()
    msfd_est = Observer(
        name="msfd_estimated",
        itoa_code="MSFD_EST",
        origin="Estimated location of CSIRO Marsfield site",
        location=est_loc,
        time=Time(60002.3, format="mjd"),
    )
    msfd_est.update()

    fitter = fit_residuals(msfd.parfile, msfd.timfile, msfd.est, msfd.err)
    fig = plot_residuals(
        fitter.resids.time_resids.to_value(u.us).astype(float),
        fitter.toas.get_mjds().value.astype(float),
        fitter.toas.get_errors().to_value(u.us).astype(float),
    )
    fig.update_layout(width=1000, height=1000)
    fig.write_image("msfd_incorrect.png")
    # fig.show()


def test_residuals_msfd_mse(msfd):
    # Vela pulse half width (m)
    half_pulsar_width = 0.08932838502359318 * SPEED_OF_LIGHT
    # establish initial true position of observers
    msfd_sim = Observer(
        name="msfd_actual",
        itoa_code="MSFD",
        location=EarthLocation.from_geocentric(
            x=-4646251.90 * u.m, y=2565112.87 * u.m, z=-3525535.83 * u.m
        ),
        origin="True location of CSIRO Marsfield telescope site",
        time=Time(60002.3, format="mjd"),
    )
    msfd_sim.update()
    msfd_est = Observer(
        name="msfd_estimated",
        itoa_code="MSFD_EST",
        origin="Estimated location of CSIRO Marsfield site",
        location=process_estimate(msfd),
        time=Time(60002.3, format="mjd"),
    )
    msfd_est.update()
    fitter = fit_residuals(msfd.parfile, msfd.timfile, msfd.est, msfd.err)
    residuals_true = fitter.resids.time_resids.to_value(u.us).astype(float)

    # Adjust position in range (x +- P0 / 2, y += P0 /2, z)
    x_range = np.linspace(msfd.est[0] - half_pulsar_width, msfd.est[0] + half_pulsar_width, 20)
    y_range = np.linspace(msfd.est[1] - half_pulsar_width, msfd.est[1] + half_pulsar_width, 20)
    residuals_mse = np.empty(shape=(len(x_range), len(y_range)))
    for i, x in enumerate(x_range):
        for j, y in enumerate(y_range):
            msfd.est = (x, y, msfd.est[2])
            est_loc = process_estimate(msfd)
            msfd_est.location = est_loc
            msfd_est.update()
            fitter = fit_residuals(msfd.parfile, msfd.timfile, msfd.est, msfd.err)
            # calculate mean squared error (MSE) of timing residuals
            mse = np.square(
                np.subtract(fitter.resids.time_resids.to_value(u.us).astype(float), residuals_true)
            ).mean()
            residuals_mse[i][j] = mse
            print(mse)

    print(residuals_mse)
    fig = plot_residuals_mse(x_range, y_range, residuals_mse)
    fig.update_layout(width=1000, height=1000)
    fig.write_image("msfd_mse_surface.png")
    # fig.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", dest="timfile", type=str, help="Path to .tim file", required=True)
    parser.add_argument("-p", dest="parfile", type=str, help="Path to .par file", required=True)
    parser.add_argument(
        "--est",
        dest="est",
        type=float,
        nargs="+",
        help="Initial position estimate for observer in WGS84",
        required=True,
    )
    parser.add_argument(
        "--err",
        dest="err",
        type=float,
        nargs="+",
        help="Initial position error for observer in WGS84",
    )
    parser.add_argument(
        "-f", dest="frame", type=str, help="Reference frame for position", required=True
    )
    args = parser.parse_args()  # noqa: F811
    test_residuals_msfd_true(args)
