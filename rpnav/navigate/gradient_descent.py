import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from numpy.typing import NDArray
from pint.fitter import Fitter, MaxiterReached
from scipy.optimize import minimize

from rpnav.observe.observer import Observer


def mse(resids: list[float]):
    return np.square(resids).mean()


def rmse(resids: list[float]):
    return np.sqrt(mse(resids))


def heuristic(
    location: tuple[float],
    observer: Observer,
    fitter: Fitter,
    maxiter: int = 100,
):
    print(location)
    observer.update(location = EarthLocation.from_geodetic(
        lon=location[0] * u.arcsec,
        lat=location[1] * u.arcsec,
    ))
    try:
        fitter.fit_toas(maxiter=maxiter)
    except MaxiterReached:
        return np.inf

    resids = mse(fitter.resids.time_resids.to_value(u.us) * 1e3)
    print(resids)
    return resids


def localise(
    position: tuple[float, float],
    observer: Observer,
    fitter: Fitter,
    tol: float = 1e-20,
):
    res = minimize(
        heuristic,
        position,
        method='BFGS',
        args=(observer, fitter),
        tol=tol
    )
    return res["x"]
