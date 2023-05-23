import astropy.units as u
import numpy as np
from astropy.coordinates import EarthLocation
from pint.fitter import Fitter, MaxiterReached
from scipy.optimize import minimize

from rpnav.observe.observer import Observer


def rmse(resids: list[float]):
    return np.sqrt(np.square(resids).mean())


def heuristic(
    location: tuple[float],
    observer: Observer,
    fitter: Fitter,
    maxiter: int = 100,
):
    observer.location = EarthLocation.from_geodetic(
        lon=location[0] * u.deg,
        lat=location[1] * u.deg,
    )
    observer.update()
    try:
        fitter.fit_toas(maxiter=maxiter)
    except MaxiterReached:
        return np.inf

    resids = fitter.resids.time_resids.to_value(u.us).astype(float)
    return rmse(resids)


def localise(
    start: tuple[float, float],
    observer: Observer,
    fitter: Fitter,
    method: str = "L-BFGS-B",
    tol: float = 1e-12,
):
    res = minimize(heuristic, start, method=method, args=(observer, fitter), tol=tol)
    return res["x"]
