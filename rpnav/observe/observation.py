from typing import Any
from math import pi, sqrt

import astropy.units as u
import numpy as np
from astropy.time import Time
from astropy.coordinates import AltAz
from astropy.coordinates.angles import Angle

from rpnav.observe import Antenna
from rpnav.observe.pulsar import Pulsar
from rpnav.conversions import jansky_to_watt
from astropy import constants as const

class Observation:
    def __init__(
        self,
        observer: Antenna,
        pulsar: Pulsar,
        localtime: Time = None,
        integration_time: u.s = None,
        snr: u.dimensionless_unscaled = None,
        receive_power: u.W = None,
        toa_err:  u.ns = None,
    ):
        self._pulsar = pulsar
        self._antenna = Antenna
        self._localtime = localtime
        self._integration_time = integration_time
        self._toa_err : u.ns = toa_err
        self._receive_power : u.W = receive_power
        self._snr : u.dimensionless_unscaled = snr
    
    def snr(self, integration_time: u.s = None) -> u.W:
        integtime = integration_time if integration_time is not None else self._integration_time
        psr = self._pulsar
        observer = self._antenna
        snr = None
        if self._snr is None:
            # snr = (Smean*jy) * Aeff * sqrt(np*deltaNu*tobs)/2.0/k/(tsys+tsky)*sqrt((psr[i].p0-w)/w);
            if self._is_set([self._receive_power, observer.noise_power]):
                snr = self._receive_power / self._noise_power
            elif self._is_set([integtime, observer.n_polarisations, observer.bandwidth, psr.pulse_peak_amplitude, 
                observer.system_temp]):
                """
                Reference:
                Buist, P.J., Engelen, S., Noroozi, A., Sundaramoorthy, P., Verhagen, S., Verhoeven, C., 
                "Principles and Potential of Pulsar Navigation," Proceedings of the 24th International 
                Technical Meeting of the Satellite Division of The Institute of Navigation (ION GNSS 2011), 
                Portland, OR, September 2011, pp. 3503-3515. 
                Equation  (15)
                """
                snr = (
                    sqrt(self._n_polarisations * integtime * observer.bandwidth) 
                    * (psr.flux_density(observer.centre_frequency) * observer.gain_dpfu / observer.system_temp)
                    * sqrt(psr.pulse_width(psr.period - psr.pulse_width))
                    / psr.period
                ) 
        else:
            snr = self._snr
        return snr
    
    def receive_power(self) -> u.W:
        psr = self._pulsar
        observer = self._antenna
        if self._is_set([observer.effective_area, observer.polarisation_param, observer.centre_frequency, observer.bandwidth]):
            receive_power = (
                observer.polarisation_param * observer._effective_area 
                * jansky_to_watt(psr.flux_density(observer.centre_frequency))
                * observer.bandwidth
            ).to(u.W)
        return receive_power
    
    def min_observable_flux_density(
        self,
        integration_time: u.s = None,
        correction: u.dimensionless_unscaled = 1,
    ) -> u.mJy:
        """Calculate minimum  narrowband spectral line sensitivity  that can be measured by the antenna
        given a fixed integration time.

        Parameters
        ----------
            integration_time (float, s):
            pulse_width_to_period (float):
            correction (float):

        Returns
        -------
            (mJ)

        Reference
        ---------
        (A. 12)
        https://link.springer.com/content/pdf/bbm:978-3-642-19627-0/1.pdf
        """
        integtime = integration_time if integration_time is not None else self._integration_time
        observer = self._antenna
        width_to_period = self._pulsar.pulse_width_10 / self._pulsar.period

        flux_dens = flux_dens
        if self._is_set([observer, observer.bandwidth, observer.radio_gain, observer.snr(integtime)]):
            flux_dens = (
                self.snr(integtime) * (correction * observer.system_temp)
                / (np.sqrt(observer.n_polarisations * observer.bandwidth * integtime) * observer.gain_dpfu)
                * np.sqrt(width_to_period / (1 - width_to_period)) 
            ).to(u.mJy)
        elif self._is_set([observer, observer.bandwidth, observer.effective_area, observer.snr(integtime)]):
            flux_dens = (
                self.snr(integtime) * (correction * observer.system_temp * 8 * const.k_B)
                / (observer.effective_area * np.sqrt(observer.n_polarisations * observer.bandwidth * integtime) * observer.gain_dpfu)
                * np.sqrt(width_to_period / (1 - width_to_period)) 
            ).to(u.mJy)
        return flux_dens

    def get_alt_az(self, localtime: Time) -> AltAz:
        """
        Get altitude and Azimuth coordinates for a pulsar given an observer location and a local time  
        
        Args:
            antenna:                Observer object
            localtime:              Time of observation

        Returns:
            Astropy AltAz object 
        """
        obstime = localtime if localtime is not None else self._localtime
        return self._to_sky_coord().transform_to(AltAz(obstime=obstime, location=self._antenna.location))
    
    def toa_err(self, integration_time: u.s = None, reference: "Observation" = None) -> u.ns:
        integtime = integration_time if integration_time is not None else self._integration_time
        observer = self._antenna
        psr = self._pulsar

        if reference.pulsar.name != psr.name:
            raise Exception("Reference pulsar must be the same as observation pulsar")
        
        if self._toa_err is None:
            if self._is_set([observer, observer.bandwidth, self.snr(integtime), reference, reference.antenna.bandwidth, reference.snr(integtime)]):
                """
                Reference:
                    Radio-Frequency Pulsar Observation using Small-Aperture Antennas (Eqn 8)
                """
                toa_err = np.sqrt(
                    reference.toa_err()**2 * (reference.snr()**2 * reference.bandwidth * reference.integration_time) 
                    / (self.snr(integtime)**2 * observer.bandwidth * integtime) 
                )
            if self._is_set([observer.bandwidth, observer.snr, reference, reference.bandwidth, reference.snr(), reference.toa_err]):
                """
                Assume observers have same integration time

                Reference:
                    Radio-Frequency Pulsar Observation using Small-Aperture Antennas (Eqn 8)
                """
                toa_err = np.sqrt((
                    (reference.toa_err**2 * reference.snr(integtime)**2 * reference.bandwidth)
                    / (self.snr(integtime)**2 * observer.bandwidth) 
                ))
            if self._is_set([self._pulse_width_10, self._period, observer.bandwidth, observer.radio_gain, observer.snr(self, integtime), integtime]):
                """
                Reference:
                    Radio-Frequency Pulsar Observation using Small-Aperture Antennas (Eqn 9)
                """
                toa_err =  (
                    np.sqrt(
                        pow(psr.pulse_width_10, 3) / (2 * psr.period)
                    ) / (
                        self.snr(integtime) * observer.radio_gain 
                        * pow(2 * pi * np.log(2), 1/4) 
                        *  np.sqrt(observer.bandwidth * integtime)
                    )
                )
        else:
            toa_err = self._toa_err
        return toa_err
    
    def is_observable(
        self,
        localtime: Time,
        integration_time: u.s,
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
        obstime = localtime if localtime is not None else self._localtime
        integtime = integration_time if integration_time is not None else self._integration_time
        observer = self._antenna
        is_vis: bool = self._pulsar.flux_density(observer.centre_frequency) > self.min_observable_flux_density(integtime)

        observable = True
        if horizon is not None:
            altaz = self.get_alt_az(observer, obstime)
            observable = altaz.alt > horizon
        return is_vis and observable
    
    @property
    def pulsar(self) -> Pulsar:
        return self._pulsar

    @property
    def antenna(self) -> Antenna:
        return self._antenna

    @property
    def integration_time(self) -> u.ms:
        return self._integration_time

    @property
    def toa_err(self) -> u.ms:
        return self._toa_err

    @property
    def localtime(self) -> Time:
        return self._localtime

    def _is_set(self, vals: list[Any]) -> bool:
        return None not in vals
    
    def _propagate_calculations(self) -> float:
        for _ in range(3):
            for var in type(self).__dict__:
                getattr(self, var)