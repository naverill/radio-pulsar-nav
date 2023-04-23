import logging
import pandas as pd
import numpy as np
import psrqpy
from pulsar_spectra.catalogue import collect_catalogue_fluxes
from pulsar_spectra.spectral_fit import find_best_spectral_fit, estimate_flux_density

from __init__ import logger
logger.disabled = True


def load_catalogue(centre_freq: float = 1400) -> pd.DataFrame:
    q = psrqpy.QueryATNF(params=[
        "NAME",     # Pulsar name
        "RAJD",     # Right ascension (J2000) (degrees) 
        "DecJD",    # Declination (J2000) (degrees)
        "GL",       # Galactic longitude (degrees) 
        "GB",       # Galactic latitude (degrees)   
        "P0",       # Barycentric period of the pulsar (s)
        "DM",       # Dispersion measure (cm-3 pc)
        "W50",      # Width of pulse at 50% of peak (ms)
        "W10",      # Width of pulse at 10% (ms)
        "S1400",    # Temp value until spectral fit is working
    ])
    #cat["S"], cat["SERR"] = np.nan, np.nan
    #cat_dict = collect_catalogue_fluxes()
    #for i, pulsar in enumerate(cat["NAME"]):
    #    if pulsar not in cat_dict:
    #        logger.warning("Pulsar config not found")
    #        continue

    #    freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
    #    if not freqs:
    #        logger.warning("Frequency values not found")
    #        continue

    #    model, m, _, _, _ = find_best_spectral_fit(
    #        pulsar, freqs, bands, fluxs, flux_errs, refs
    #    )
    #    if model:
    #        cat["S"].iloc[i], cat["SERR"].iloc[i] = estimate_flux_density(centre_freq, model, m)
    #    else:
    #        logger.warning(f"Failed to find best fit spectral model for {pulsar}")
    return q.table

