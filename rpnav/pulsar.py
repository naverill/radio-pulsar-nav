import pandas as pd
import psrqpy
from astropy.coordinates import ICRS, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time

from rpnav import logger
from rpnav.antenna import Antenna

# from pulsar_spectra.catalogue import collect_catalogue_fluxes
# from pulsar_spectra.spectral_fit import estimate_flux_density, find_best_spectral_fit


logger.disabled = True


class Pulsar(SkyCoord):
    def __init__(
        self,
        name: str,
        centre_freq: float,
        flux_density: float,
        ra: float,
        dec: float,
        period: float,
        dispersion_measure: float,
        pulse_width_50: float,
        pulse_width_10: float,
        flux_density_err: float = None,
    ):
        """
        "NAME",  # Pulsar name
        "RAJD",  # Right ascension (J2000) (degrees)
        "DecJD",  # Declination (J2000) (degrees)
        "P0",  # Barycentric period of the pulsar (s)
        "DM",  # Dispersion measure (cm-3 pc)
        "W50",  # Width of pulse at 50% of peak (ms)
        "W10",  # Width of pulse at 10% (ms)
        """
        self.name = name
        self.flux_density = flux_density
        self.period = period
        self.dispersion_measure = dispersion_measure
        self.pulse_width_50 = pulse_width_50
        self.pulse_width_10 = pulse_width_10
        super().__init__(frame=ICRS, ra=ra, dec=dec, unit="deg")

    def to_sky_coord(self) -> SkyCoord:
        return SkyCoord(frame=ICRS, ra=self.ra, dec=self.dec, unit="deg")

    def is_observable(
        self,
        antenna: Antenna,
        t: Time,
        integ_time: float,
        horizon: float = 0,
    ) -> bool:
        is_vis: bool = self.flux_density > antenna.min_observable_flux_density(integ_time)
        observable = True
        if horizon is not None:
            altaz = self.get_alt_az(antenna, t)
            observable = altaz.alt > horizon
        return is_vis and observable

    def get_alt_az(self, antenna: Antenna, t: Time):
        return self.to_sky_coord().transform_to(AltAz(obstime=t, location=antenna.location))

    @staticmethod
    def load_catalogue(centre_freq: float = 1400) -> pd.DataFrame:
        """Load pulsar data from PSR Catalogue."""
        cat = psrqpy.QueryATNF(
            params=[
                "NAME",  # Pulsar name
                "RAJD",  # Right ascension (J2000) (degrees)
                "DECJD",  # Declination (J2000) (degrees)
                "P0",  # Barycentric period of the pulsar (s)
                "DM",  # Dispersion measure (cm-3 pc)
                "W50",  # Width of pulse at 50% of peak (ms)
                "W10",  # Width of pulse at 10% (ms)
                "S1400",
            ]
        ).table
        pulsars: list[Pulsar] = []
        # cat_dict = collect_catalogue_fluxes()
        for row in cat:
            flux_dens = row["S1400"]
            flux_dens_err = None
            #    if pid not in cat_dict:
            #        logger.warning("Pulsar config not found")
            #        continue

            #    freqs, bands, fluxs, flux_errs, refs = cat_dict[pid]
            #    if not freqs:
            #        logger.warning("Frequency values not found")
            #        continue
            #    model, m, _, _, _ = find_best_spectral_fit(pulsar, freqs, bands, fluxs, flux_errs, refs)
            #    if model:
            #        flux_dens, flux_dens_err = estimate_flux_density(
            #            centre_freq, model, m
            #        )
            #    else:
            #        logger.warning(f"Failed to find best fit spectral model for {pulsar}")
            if not flux_dens:
                continue
            pulsars.append(
                Pulsar(
                    name=row["NAME"],
                    ra=row["RAJD"],
                    dec=row["DECJD"],
                    period=row["P0"],
                    dispersion_measure=row["DM"],
                    pulse_width_50=row["W50"],
                    pulse_width_10=row["W10"],
                    flux_density=flux_dens,
                    flux_density_err=flux_dens_err,
                    centre_freq=centre_freq,
                )
            )
        return pulsars
