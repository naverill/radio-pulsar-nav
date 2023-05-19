from typing import TextIO, Union

import astropy.units as u
from pint.fitter import Fitter
from pint.models import get_model, get_model_and_toas


def fit_residuals(
    parfile: Union[TextIO, str], timfile: Union[TextIO, str], pos: list[float], err: list[float]
):
    model, toas = get_model_and_toas(parfile, timfile)
    fitter = Fitter.auto(toas, model)
    fitter.fit_toas(maxiter=100)
    return fitter
