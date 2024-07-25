from enum import Enum, auto
from typing import TextIO, Union

import astropy.units as u
from pint import logging
from pint.fitter import Fitter, MaxiterReached
from pint.models import get_model, get_model_and_toas


def fit_residuals(
    parfile: Union[TextIO, str],
    timfile: Union[TextIO, str],
    fitter: Fitter = None,
    maxiter: int = 100,
):
    """
    Fit timing model based on input parameter and timing file

    Args:
        parfile:    Pulsar parameter file (.par)
        timfile:    Pulsar timing observation file (.tim)
        fitter:     Object to fit timing models to TOAs can be 
                    pint.fitter.WLSFitter for basic fitting, 
                    pint.fitter.GLSFitter for fitting with noise
                    models that imply correlated errors, and 
                    pint.fitter.WidebandTOAFitter for TOAs that 
                    contain DM information. If none is supplied,
                    automatic fitting is applied 
        maxiter:    Maximum number of steps allowed in fitting 
                    process

    Returns:
        Fitted timing model
    """
    model, toas = get_model_and_toas(parfile, timfile)
    if fitter is None:
        fit = Fitter.auto(toas, model)
    else:
        fit = fitter(toas=toas, model=model)

    try:
        fit.fit_toas(maxiter=maxiter)
    except MaxiterReached:
        logging.error("PINT Residuals Fitter failed to fully converge.")
    return fit
