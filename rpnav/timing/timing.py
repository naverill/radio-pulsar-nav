from enum import Enum, auto
from typing import TextIO, Union

import astropy.units as u
from pint import logging
from pint.fitter import Fitter, MaxiterReached
from pint.models import get_model, get_model_and_toas
from loguru import logger as log

logging.setup()


def fit_residuals(
    parfile: Union[TextIO, str],
    timfile: Union[TextIO, str],
    fitter: Fitter = None,
    maxiter: int = 100,
) -> Union[Fitter, None]:
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
    model, toas = get_model_and_toas(parfile, timfile, allow_T2=True, allow_tcb=True)
    if fitter is None:
        fit = Fitter.auto(toas, model)
    else:
        fit = fitter(toas=toas, model=model)

    try:
        fit.fit_toas(maxiter=maxiter)
    except MaxiterReached:
        log.warning("PINT Residuals Fitter failed to fully converge.")
        return None
    return fit
