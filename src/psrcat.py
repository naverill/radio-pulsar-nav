import pandas as pd
import psrqpy
from pulsar_spectra.catalogue import collect_catalogue_fluxes
from pulsar_spectra.spectral_fit import find_best_spectral_fit


def load_catalogue(centre_freq: float) -> pd.DataFrame:
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
    cat = q.pandas
    #cat_dict = collect_catalogue_fluxes()
    #for i, pulsar in zip(cat.index, cat["NAME"]):
    #    freqs, bands, fluxs, flux_errs, refs = cat_dict[pulsar]
    #    model, m, _, _, _ = find_best_spectral_fit(
    #        pulsar, freqs, bands, fluxs, flux_errs, refs
    #    )
    #    cat[i, "S"], cat[i, "SERR"] = estimate_flux_density(centre_freq, model, m)
    return cat

