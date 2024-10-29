import pandas as pd
import psrqpy
from typing import Any

import astropy.units as u
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
        ref_frequency: u.MHz,
        ref_flux_density: u.mJy,
        ra: Angle,
        dec: Angle,
        period: u.ms,
        distance: u.kpc,
        dispersion_measure: u.cm**-3 * u.pc,
        pulse_width_50: u.ms,
        pulse_width_10: u.ms,
        spectral_index: u.dimensionless_unscaled = None,
    ):
        """
        Initialise Pulsar object

        Parameters
        ----------
            namw:                   Pulsar name
            RAJD:                   Right ascension (J2000) (degrees)
            DecJD:                  Declination (J2000) (degrees)
            period:                 Barycentric period of the pulsar (s)
            dispersion_measure:     Dispersion measure (cm-3 pc)
            pulse_width_50:         Width of pulse at 50% of peak (ms)
            pulse_width_10:         Width of pulse at 10% (ms)
            ref_frequency:          reference frequence with extensive measurement (MHz)
            ref_flux_density:       reference flux density - brightness at known measurement frequency (mJy)
        """
        self._name = name
        self._period = period
        self._dispersion_measure = dispersion_measure
        self._pulse_width_50 = pulse_width_50
        self._pulse_width_10 = pulse_width_10
        self._ref_frequency = ref_frequency
        self._distance = distance
        self._ref_flux_density = ref_flux_density
        self._spectral_index = spectral_index
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
        integtime: u.s,
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
        is_vis: bool = self.flux_density(antenna.centre_frequency) > antenna.min_observable_flux_density(integtime)

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
        return self._to_sky_coord().transform_to(AltAz(obstime=localtime, location=antenna.location))

    def interpolate_flux_density(
        self, centre_freq: u.MHz
    ) -> tuple[u.mJy, u.mJy]:
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
        cat_dict = collect_catalogue_fluxes()

        params = [
            "NAME",
            "RAJD",             # Right ascension (J2000) (degrees)
            "DECJD",            # Declination (J2000) (degrees)
            "P0",               # Barycentric period of the pulsar (s)
            "DM",               # Barycentric period of the pulsar (s)
            "W50",              # Width of pulse at 50% of peak (ms).
            "W10",              # Width of pulse at 10% (ms).
            "DIST",             # est estimate of the pulsar distance using the YMW16 DM-based distance as default (kpc).
            "SPINDX",           # Radio spectral index
        ]
        # Query PSRCAT to get complete list of pulsar profiles
        try:
            cat = psrqpy.QueryATNF(params=params, psrs=[self.name]).table
        except ValueError:
            raise Exception(
                "Failed to read from psrcat. Check that centre frequency" " is a valid table option"
            )
        
        for row in cat:
            if row["NAME"] != self.name:
                continue
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
                    flux_dens, flux_dens_err = estimate_flux_density(centre_freq.to_value(u.MHz), model, m)
                else:
                    logger.warning(f"Failed to find best fit spectral model for {row['NAME']}")

            if flux_dens is not None:
                flux_dens = flux_dens * u.mJy

            if flux_dens_err is not None:
                flux_dens = flux_dens_err * u.mJy

        return flux_dens, flux_dens_err

    @staticmethod
    def load_catalogue(ref_freq: u.MHz = 1400 * u.MHz) -> list["Pulsar"]:
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
        valid_freq = [400,1400,2000,30,40,50,60,80,100,150,200,300,350,600,700,800,900,1600,3000,4000,6000,8000]

        if ref_freq is None or ref_freq.to_value(u.MHz) not in valid_freq:
            print(f"Invalid reference frequency. must be in [{', '.join(str(i) for i in valid_freq)}]. Defaulting to 1400 MHz")
            ref_freq = 1400 * u.MHz

        scode = f"S{int(ref_freq.to_value(u.MHz))}"

        params = [
            "NAME",
            "RAJD",             # Right ascension (J2000) (degrees)
            "DECJD",            # Declination (J2000) (degrees)
            "P0",               # Barycentric period of the pulsar (s)
            "DM",               # Barycentric period of the pulsar (s)
            "W50",              # Width of pulse at 50% of peak (ms).
            "W10",              # Width of pulse at 10% (ms).
            "DIST",             # est estimate of the pulsar distance using the YMW16 DM-based distance as default (kpc).
            "SPINDX",           # Radio spectral index
            scode,              # Mean flux density at reference frequency MHz (mJy)
        ]

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

            # Create pulsar object from profile
            pulsars.append(
                Pulsar(
                    name=row.get("NAME"),
                    ra=row.get("RAJD") * u.deg,
                    dec=row.get("DECJD") * u.deg,
                    period=row.get("P0") * u.s,
                    dispersion_measure=row.get("DM") * u.cm**-3 * u.pc,
                    pulse_width_50=row.get("W50") * u.ms,
                    pulse_width_10=row.get("W10") * u.ms,
                    # flux_density_err=flux_dens_err * u.mJy if flux_dens_err is not None else None,
                    ref_frequency=ref_freq,
                    ref_flux_density=row.get(scode) * u.mJy,
                    spectral_index=row.get("SPINDX"),
                    distance=row.get("DIST"),
                )
            )

        return pulsars

    @property
    def timing_uncertainty(self) -> u.ms:
        """
        Estimate for intrinsic timing error for the pulsar. This is calculated as one tenth of the
        Width of pulse at 10%
        
        Args:
            None

        Returns:
            Puldar timing uncertainty
        """
        return self.pulse_width_10 / 10

    @property
    def name(self) -> str:
        return self._name 

    def flux_density(self, centre_frequency: u.MHz) -> u.mJy:
        """
        Flux density of the pulsar at centre frequency

        Reference:
            https://www.researchgate.net/publication/280085456_Radio_Pulsar_Receiver_Systems_for_Space_Navigation
        """
        if self._is_set([self._ref_frequency, self._spectral_index, self._ref_flux_density]):
            flux_density = (self._ref_flux_density 
                * (centre_frequency.to(u.GHz) / self._ref_frequency.to(u.GHz))**(self._spectral_index)
            )
        return flux_density 

    @property
    def flux_density_err(self) -> u.mJy:
        return self._flux_density_err

    @property
    def dispersion_measure(self) -> u.MHz:
        return self._dispersion_measure

    @property
    def period(self) -> u.ms:
        return self._period

    @property
    def pulse_width_50(self) -> u.ms:
        return self._pulse_width_50

    @property
    def pulse_width_10(self) -> u.ms:
        return self._pulse_width_10

    @property
    def distance(self) -> u.kpc:
        return self._distance

    @property
    def spectral_index(self) -> u.dimensionless_unscaled:
        return self._spectral_index
    
    def _is_set(self, vals: list[Any]) -> bool:
        return None not in vals
    
    def _propagate_calculations(self) -> float:
        for _ in range(3):
            for var in type(self).__dict__:
                getattr(self, var)
