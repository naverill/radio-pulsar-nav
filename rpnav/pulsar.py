import pandas as pd
import psrqpy
from astropy.coordinates import ICRS, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time
from astropy.coordinates.angles import Angle
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
        Initialise Pulsar object

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
        """
        Convert Pulsar object to Astropy SkyCoord object
        
        Args:
            None

        Returns:
            SkyCoord object in the ICRS frame 
        """
        return SkyCoord(frame=ICRS, ra=self.ra, dec=self.dec, unit="deg")

    def is_observable(
        self,
        antenna: Antenna,
        localtime: Time,
        integtime: float,
        horizon: Angle = None,
    ) -> bool:
        """
        Calculate whether pulsar is observable by an observer at a specific time. To be 
        observable the pulsar must have a flux density above a certain threshold for the 
        antenna configuration. If a horizon value is submitted, the function will 
        additinally return if the observer has line-of-sight to the pulsar. 
        
        Args:
            antenna:                Observer object
            t:                      Time of observation
            integ_time:             Length of observation time
            horizon (optional):     Minumum altitude angle to be considered visible


        Returns:
            Boolean indicating if pulsar is visible 
        """
        is_vis: bool = self.flux_density > antenna.min_observable_flux_density(integtime)

        observable = True
        if horizon is not None:
            altaz = self.get_alt_az(antenna, localtime)
            observable = altaz.alt > horizon
        return is_vis and observable

    def get_alt_az(self, antenna: Antenna, localtime: Time) -> AltAz:
        """
        Get altitude and Azimuth coordinates for a pulsar given an observer location and a local time  
        
        Args:
            antenna:                Observer object
            localtime:              Time of observation

        Returns:
            Astropy AltAz object 
        """
        return self._to_sky_coord().transform_to(AltAz(obstime=t, location=antenna.location))

    @staticmethod
    def _interpolate_flux_density(
        centre_freq: int, row: Table, cat_dict: dict
    ) -> tuple[float, float]:
        """
        Estimate the flux density for a pulsar at a specific centre frequency. Pulsars have been observed 
        at limited bandwidths, and the flux density of the pulsars do not vary linearly with the centre frequency [1].
        Therefore, to estimate the radio luminosity at an unobserved frequency, we must interpolate given previous observations
        using the pulsar_spectra package [2]. This package is released and maintained by Swainston et al. 
        
        Args:
            centre_freq:            Centre frequency of observation
            row:                    Pulsar profile information
            cat_dict:               Catalogue data containing all recorded flux density observations for a i=given pulsar

        Returns:
            Tuple containing interpolared flux density at centre frequency, and estimated error  

        Citations
            [1] Swainston, N. A., C. P. Lee, S. J. McSweeney, and N. D. R. Bhat. “Pulsar_spectra: A Pulsar Flux Density Catalogue and Spectrum Fitting Repository.” 
            Publications of the Astronomical Society of Australia 39 (2022): e056. https://doi.org/10.1017/pasa.2022.52. 
            [2] pulsar_spectra Github, https://github.com/NickSwainston/pulsar_spectra
        """
        flux_dens = None
        flux_dens_err = None

        # Pulsar 
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
        """
        Load pulsar data from PSR Catalogue (PSRCAT).
        
        Args:
            centre_freq:    centre frequency of observation in MHz
            interpolate:    Boolean indicating if pulsar flux density should be interpolated for 
                            the given centre frequency. If false, the catalogue will check directly
                            for flux density values at the given centre frequency, and fail if that
                            frequency has not been observed 

        Returns:
            List of Pulsar objects for every pulsar in ANTF catalogue
        """
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

        # Query PSRCAT to get complete list of pulsar profiles
        try:
            cat = psrqpy.QueryATNF(params=params).table
        except ValueError:
            raise Exception(
                "Failed to read from psrcat. Check that centre frequency" " is a valid table option"
            )

        # Extract pulsar profile for selected centre frequency 
        pulsars: list[Pulsar] = []
        for row in cat:
            flux_dens = None
            flux_dens_err = None

            # Interpolate centre frequency if requested
            if interpolate:
                flux_dens, flux_dens_err = Pulsar._interpolate_flux_density(
                    centre_freq, row, cat_dict
                )
            # Check for specific centre frequency
            else:
                flux_dens = row.get(f"S{centre_freq}")

            if not flux_dens:
                continue

            # Create pulsar object from profile
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
    def timing_uncertainty(self) -> float:
        """
        Estimate for intrinsic timing error for the pulsar. This is calculated as one tenth of the
        Width of pulse at 10%
        
        Args:
            None

        Returns:
            Puldar timing uncertainty
        """
        return self.pulse_width_10 / 10
