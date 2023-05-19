import argparse
import io
import json
from math import pi
from pathlib import Path

import astropy.units as u
import pytest
from astropy.coordinates import EarthLocation
from astropy.time import Time
from pint.observatory.topo_obs import load_observatories

from rpnav.observe.antenna import Antenna
from rpnav.observe.observer import Observer
from rpnav.timing.timing import fit_residuals
from rpnav.timing.visualise import plot_residuals

FILE_DIR = Path(__file__).parent


class Args:
    sim = "parkes_sim"
    timfile = f"{FILE_DIR}/../simulate/parkes_sim/output/real_0/J0835-4510.tim"
    parfile = f"{FILE_DIR}/../simulate/parkes_sim/output/real_0/J0835-4510.tdb.par"
    est = [-4460892.6, 3682358.9, -3696756.0]
    frame = "geocentric"
    err = None


@pytest.fixture(scope="module")
def args():
    return Args()


def test_residuals_parkes(args):
    est_loc = process_estimate(args)
    pks_sim = Observer(
        name="pks",
        itoa_code="PARKES",
        location=est_loc,
        origin="Actual location of Parkes telescope site",
        time=Time(60002.3, format="mjd"),
    )
    load_observatories(io.StringIO(json.dumps(pks_sim.to_json())), overwrite=True)
    fitter = fit_residuals(args.parfile, args.timfile, args.est, args.err)
    fig = plot_residuals(
        fitter.resids.time_resids.to_value(u.us).astype(float),
        fitter.toas.get_mjds().value.astype(float),
        fitter.toas.get_errors().to_value(u.us).astype(float),
    )
    fig.update_layout(width=1000, height=1000)
    fig.write_image("parkes.png")
    # fig.show()


def test_residuals_msfd(args):
    est_loc = process_estimate(args)
    msfd_sim = Observer(
        name="msfd_actual",
        itoa_code="MSFD",
        #        location=est_loc,
        location=EarthLocation.from_geocentric(
            x=-4460892.6 * u.m, u=3682358.9 * u.m, z=-3696756.0 * u.m
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
    fitter = fit_residuals(args.parfile, args.timfile, args.est, args.err)
    fig = plot_residuals(
        fitter.resids.time_resids.to_value(u.us).astype(float),
        fitter.toas.get_mjds().value.astype(float),
        fitter.toas.get_errors().to_value(u.us).astype(float),
    )
    fig.update_layout(width=1000, height=1000)
    fig.write_image("test.png")
    # fig.show()


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
