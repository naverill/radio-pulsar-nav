import psrqpy
from pulsar_spectra.catalogue import collect_catalogue_fluxes
from pulsar_spectra.spectral_fit import find_best_spectral_fit


def load_catalogue():
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
    ])

    # get frequencies as an astropy table
    return q.pandas
