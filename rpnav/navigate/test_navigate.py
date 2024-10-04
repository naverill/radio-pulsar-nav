from pathlib import Path
import csv

import pandas as pd
import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation
from astropy.time import Time
from pint import logging
from pint.fitter import DownhillWLSFitter
from pint.models import get_model_and_toas

from rpnav.constants import SPEED_OF_LIGHT
from rpnav.navigate.gradient_descent import localise, rmse
from rpnav.navigate.visualise import plot_heatmap
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
    sim = "msfd_est_sim"
    loc = EarthLocation.from_geocentric(
        #        lat=-33.77290643916046 * u.deg, lon=151.0976937264337 * u.deg
        x=-4646251.90 * u.m,
        y=2565112.87 * u.m,
        z=-3525535.83 * u.m,
    ).geocentric
    return Args(
        timfile=f"{FILE_DIR}/../simulate/{sim}/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/../simulate/{sim}/output/real_0/J0835-4510.tdb.par",
        est=[loc[0].to_value(u.m), loc[1].to_value(u.m), loc[2].to_value(u.m)],
        frame="geocentric",
        err=None,
    )

def test_gradient_descent(msfd):
    half_pulse_width = (0.08932838502359318 * u.s / 2) * SPEED_OF_LIGHT
    loc_true = process_estimate(msfd)
    N_ITER = 100
    msfd_sim = Observer(
        name="msfd_actual",
        itoa_code="MSFD",
        location=loc_true,
        origin="True location of CSIRO Marsfield telescope site",
        time=Time(60002.3, format="mjd"),
    )
    msfd_est = Observer(
        name="msfd_estimated",
        itoa_code="MSFD_EST",
        origin="Estimated location of CSIRO Marsfield site",
        location=loc_true,
        time=Time(60002.3, format="mjd"),
    )
    model, toas = get_model_and_toas(msfd.parfile, msfd.timfile)
    fitter = DownhillWLSFitter(toas, model)
    fitter.fit_toas()
    resids = fitter.resids.time_resids.to_value(u.us)
    initial_err = np.empty(shape=(N_ITER, 2))
    fitted_err = np.empty(shape=(N_ITER, 2))
    for i in range(N_ITER):
        print(i)
        loc_err = EarthLocation.from_geocentric(
            x=np.random.uniform(
                low=(loc_true.x - half_pulse_width).to_value(u.m),
                high=(loc_true.x + half_pulse_width).to_value(u.m),
            )
            * u.m,
            y=np.random.uniform(
                low=(loc_true.y - half_pulse_width).to_value(u.m),
                high=(loc_true.y + half_pulse_width).to_value(u.m),
            )
            * u.m,
            z=loc_true.z,
        )
        msfd_est.update(location=loc_err)
        start = (loc_err.geodetic.lon.to_value(u.arcsec), loc_err.geodetic.lat.to_value(u.arcsec))
        loc_est = localise(start, msfd_est, fitter)
        loc_est = EarthLocation.from_geodetic(lon=loc_est[0] * u.arcsec, lat=loc_est[1] * u.arcsec)
        initial_err[i][0] = (loc_err.lon - loc_true.lon).to_value(u.arcmin)
        initial_err[i][1] = (loc_err.lat - loc_true.lat).to_value(u.arcmin)
        fitted_err[i][0] = (loc_est.lon - loc_true.lon).to_value(u.arcmin)
        fitted_err[i][1] = (loc_est.lat - loc_true.lat).to_value(u.arcmin)
    print("Fitted errors", np.mean(fitted_errors, axis=0))


def test_rmse():
    fpath = "../simulate/outputs/navigate_sim/output/real_0/"

    df1 = pd.read_csv(fpath + "ALPSMLC30_S033E147_DSM_rmse.csv")
    df1 = pd.read_csv(fpath + "ALPSMLC30_S033E148_DSM_rmse.csv")
    df1 = pd.read_csv(fpath + "ALPSMLC30_S033E149_DSM_rmse.csv")  
    df1.set_index(np.linspace(33, 34, 3599, endpoint=False), inplace=True)
    df1.columns = np.linspace(148, 149, 3601, endpoint=False)


    df2 = pd.read_csv(fpath + "ALPSMLC30_S034E147_DSM_rmse.csv")
    df2 = pd.read_csv(fpath + "ALPSMLC30_S034E148_DSM_rmse.csv")
    df2 = pd.read_csv(fpath + "ALPSMLC30_S034E149_DSM_rmse.csv")
    df2.set_index(np.linspace(34, 35, 3599, endpoint=False), inplace=True)
    df2.columns = np.linspace(148, 150, 3601*2, endpoint=False) 

    df1 = pd.concat([df1, df2], axis=0, join="outer")
    fig = plot_heatmap(df1)

    # fig = go.Figure(data=go.Heatmap(z=rmseMap),
    #     color_scale='RdBu_r', origin='lower')
    fig.write_image("outputs/residErr8.png")