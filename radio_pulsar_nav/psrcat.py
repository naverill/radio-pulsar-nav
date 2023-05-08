"""
PSRCAT interface module
"""
import pandas as pd
import psrqpy
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.time import Time
from pulsar_spectra.catalogue import collect_catalogue_fluxes
from pulsar_spectra.spectral_fit import estimate_flux_density, find_best_spectral_fit

from radio_pulsar_nav import logger
from radio_pulsar_nav.antenna import Antenna

logger.disabled = True


def load_catalogue(centre_freq: float = 1400) -> pd.DataFrame:
    """Load pulsar data from PSR Catalogue."""
    cat = psrqpy.QueryATNF(
        params=[
            "NAME",  # Pulsar name
            "RAJD",  # Right ascension (J2000) (degrees)
            "DecJD",  # Declination (J2000) (degrees)
            "GL",  # Galactic longitude (degrees)
            "GB",  # Galactic latitude (degrees)
            "P0",  # Barycentric period of the pulsar (s)
            "DM",  # Dispersion measure (cm-3 pc)
            "W50",  # Width of pulse at 50% of peak (ms)
            "W10",  # Width of pulse at 10% (ms)
            "S1400",
        ]
    )
    flux_cat = pd.DataFrame(index=cat["NAME"], columns=["S", "SERR"])
    cat_dict = collect_catalogue_fluxes()
    for pulsar in cat["NAME"]:
        if pulsar not in cat_dict:
            logger.warning("Pulsar config not found")
            continue

        freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
        if not freqs:
            logger.warning("Frequency values not found")
            continue
        flux_cat["S"] = cat["S1400"]
        # model, m, _, _, _ = find_best_spectral_fit(pulsar, freqs, bands, fluxs, flux_errs, refs)
        # if model:
        #    flux_cat["S"][pulsar], flux_cat["SERR"][pulsar] = estimate_flux_density(
        #        centre_freq, model, m
        #    )
        # else:
        #    logger.warning(f"Failed to find best fit spectral model for {pulsar}")
    return cat.table, flux_cat


def calculate_observable_pulsars(
    cat: Table,
    flux_cat: pd.DataFrame,
    antenna: Antenna,
    t: Time,
    integ_time: float,
    above_horiz: bool = True,
) -> list[bool]:
    vis: list[bool] = []
    for i in range(len(cat)):
        sky_pos = SkyCoord(frame="galactic", l=cat["GL"][i], b=cat["GB"][i], unit="deg")
        is_vis: bool = flux_cat["S"][i] > antenna.min_observable_flux_density(integ_time)
        if above_horiz:
            is_vis = antenna.target_is_up(t, sky_pos) and is_vis
        vis.append(is_vis)
    return vis
