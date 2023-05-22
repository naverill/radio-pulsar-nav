import pandas as pd
import psrqpy
from astropy.coordinates import ICRS, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time
from pulsar_spectra.catalogue import collect_catalogue_fluxes
from pulsar_spectra.spectral_fit import estimate_flux_density, find_best_spectral_fit

from rpnav import logger
from rpnav.observe.antenna import Antenna

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
        Parameters
        ----------
            NAME:   Pulsar name
            RAJD:   Right ascension (J2000) (degrees)
            DecJD:  Declination (J2000) (degrees)
            P0:     Barycentric period of the pulsar (s)
            DM:     Dispersion measure (cm-3 pc)
            W50:    Width of pulse at 50% of peak (ms)
            W10:    Width of pulse at 10% (ms)
        """
        self.name = name
        self.flux_density = flux_density
        self.period = period
        self.dispersion_measure = dispersion_measure
        self.pulse_width_50 = pulse_width_50
        self.pulse_width_10 = pulse_width_10
        super().__init__(frame=ICRS, ra=ra, dec=dec, unit="deg")

    def _to_sky_coord(self) -> SkyCoord:
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
        return self._to_sky_coord().transform_to(AltAz(obstime=t, location=antenna.location))

    @staticmethod
    def _interpolate_flux_density(
        centre_freq: int, row: Table, cat_dict: dict
    ) -> tuple[float, float]:
        flux_dens = None
        flux_dens_err = None
        if row["NAME"] not in cat_dict:
            logger.warning("Pulsar config not found")
        else:
            model = None
            freqs, bands, fluxs, flux_errs, refs = cat_dict[row["NAME"]]
            if not freqs:
                logger.warning("Frequency values not found")
            else:
                model, m, _, _, _ = find_best_spectral_fit(
                    row["NAME"], freqs, bands, fluxs, flux_errs, refs
                )
            if model:
                flux_dens, flux_dens_err = estimate_flux_density(centre_freq, model, m)
            else:
                logger.warning(f"Failed to find best fit spectral model for {row['NAME']}")
        return flux_dens, flux_dens_err

    @staticmethod
    def load_catalogue(centre_freq: int = 1400, interpolate=False) -> list["Pulsar"]:
        """Load pulsar data from PSR Catalogue (PSRCAT)."""
        cat_dict = None
        params = [
            "NAME",
            "RAJD",
            "DECJD",
            "P0",
            "DM",
            "W50",
            "W10",
        ]
        if interpolate:
            cat_dict = collect_catalogue_fluxes()
        else:
            params.append(f"S{centre_freq}"),

        try:
            cat = psrqpy.QueryATNF(params=params).table
        except ValueError:
            raise Exception(
                "Failed to read from psrcat. Check that centre frequency" " is a valid table option"
            )

        pulsars: list[Pulsar] = []
        for row in cat:
            flux_dens = None
            flux_dens_err = None

            if interpolate:
                flux_dens, flux_dens_err = Pulsar._interpolate_flux_density(
                    centre_freq, row, cat_dict
                )
            else:
                flux_dens = row[f"S{centre_freq}"]

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

    @property
    def timing_uncertainty(self):
        return self.pulse_width_10 / 10
