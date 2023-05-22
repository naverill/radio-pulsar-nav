import argparse
import io
import json
from math import pi
from pathlib import Path

import astropy.units as u
import pytest
from astropy.coordinates import EarthLocation
from astropy.time import Time
from pint.observatory import get_observatory
from pint.observatory.topo_obs import load_observatories

from rpnav.observe.antenna import Antenna
from rpnav.observe.observer import Observer
from rpnav.timing.timing import fit_residuals
from rpnav.timing.visualise import plot_residuals

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


def test_residuals_parkes(parkes):
    est_loc = process_estimate(parkes)
    pks_sim = Observer(
        name="pks",
        itoa_code="PARKES",
        location=est_loc,
        origin="Actual location of Parkes telescope site",
        time=Time(60002.3, format="mjd"),
    )
    load_observatories(io.StringIO(json.dumps(pks_sim.to_json())), overwrite=True)
    fitter = fit_residuals(parkes.parfile, parkes.timfile, parkes.est, parkes.err)
    fig = plot_residuals(
        fitter.resids.time_resids.to_value(u.us).astype(float),
        fitter.toas.get_mjds().value.astype(float),
        fitter.toas.get_errors().to_value(u.us).astype(float),
    )
    fig.update_layout(width=1000, height=1000)
    fig.write_image("parkes.png")
    # fig.show()


def test_residuals_msfd(msfd):
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
    msfd_est = Observer(
        name="msfd_estimated",
        itoa_code="MSFD_EST",
        origin="Estimated location of CSIRO Marsfield site",
        location=est_loc,
        time=Time(60002.3, format="mjd"),
    )
    obs = msfd_sim.to_json()
    obs.update(msfd_est.to_json())

    load_observatories(io.StringIO(json.dumps(obs)), overwrite=True)
    fitter = fit_residuals(msfd.parfile, msfd.timfile, msfd.est, msfd.err)
    fig = plot_residuals(
        fitter.resids.time_resids.to_value(u.us).astype(float),
        fitter.toas.get_mjds().value.astype(float),
        fitter.toas.get_errors().to_value(u.us).astype(float),
    )
    fig.update_layout(width=1000, height=1000)
    fig.write_image("test.png")
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
    test_residuals_msfd(args)
