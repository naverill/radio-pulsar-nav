"""
Antenna observer object
"""
from math import pi, sqrt
from typing import Any

import astropy.units as u
from astropy import constants as const

import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

from rpnav.constants import BACKGROUND_TEMP
from rpnav.conversions import frequency_to_wavelength, wavelength_to_frequency, joule_to_jansky, jansky_to_watt
from rpnav.observe.observer import Observer


class Antenna(Observer):
    """
    Antenna observer object
    """

    def __init__(
        self,
        name: str,
        system_temp:  u.K,
        bandwidth: u.MHz,
        location: EarthLocation,
        snr: u.dimensionless_unscaled = None,
        aliases: list[str] = [],
        timescale: str = "TCB",
        n_polarisations: int = 1,
        solar_distance: u.AU = 1,
        time: Time = None,
        itoa_code: str = None,
        centre_frequency: u.MHz = None,
        wavelength: u.m = None,
        effective_area: u.m * u.m = None,
        diameter: u.m = None,
        n_element_coherent: int = None,
        n_element_incoherent: int = None,
        gain_dpfu: u.K / u.Jy = None,
        astronomy_gain: u.m * u.m / u.Jy = None,
        elevation: u.m = None,
        noise_power: u.W = None,
        solar_noise_temp:  u.K = None,
        galaxy_noise_temp:  u.K = None,
        side_lobe_attenuation:  u.dB = None,
        toa_err:  u.ns = None,
    ):
        """Instantiate antenna.

        Parameters
        ----------
            signal_to_noise (float): Signal to noise ratio
            system_temp (float, K): System temperature
            gain (float, K/Jy): Antenna gain
            bandwidth (float, Hz): Antenna bandwidth
            centre_freq (float, MHz):
            half_beamwidth (float, degrees):
            diameter: (float, m)
            n_element_coherent (float):
            n_element_incoherent (float):
        """
        self._snr : u.dimensionless_unscaled = snr
        self._system_temp : u.K = system_temp
        self._bandwidth : u.MHz = bandwidth
        self._centre_frequency : u.MHz = centre_frequency
        self._wavelength : u.m = wavelength
        self._gain_dpfu = gain_dpfu
        self._effective_area : u.m * u.m = effective_area
        self._n_element_coherent = n_element_coherent
        self._n_element_incoherent = n_element_incoherent
        self._diameter : u.m = diameter
        self._noise: u.W = noise_power
        self._n_polarisations: int = n_polarisations
        self._elevation: u.m = elevation
        self._location : EarthLocation = location
        self._solar_distance = solar_distance
        self._solar_noise_temp = solar_noise_temp
        self._galaxy_noise_temp = galaxy_noise_temp
        self._noise_power = noise_power
        self._side_lobe_attenuation = side_lobe_attenuation
        self._gain_dpfu_eqn = ""
        self._radio_gain_eqn = ""
        self._radio_gain = None
        self._polsarisation_param = None
        self._pulse_peak_amplitude = None
        self._toa_err : u.ns = toa_err
        self._propagate_calculations()
        super().__init__(name, location, aliases=aliases, timescale=timescale, time=time, origin="test", itoa_code=name.upper())

    def min_observable_flux_density(
        self,
        integration_time: u.s,
        pulse_width_to_period: u.dimensionless_unscaled = 0.1,
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
        return (
            self.snr() * (correction * self.system_temp)
            / (np.sqrt(self.n_polarisations * self.bandwidth * integration_time) * self.gain_dpfu)
            * np.sqrt(pulse_width_to_period / (1 - pulse_width_to_period)) 
        ).to(u.mJy)
    
    @property
    def gain_dpfu(self) -> u.K / u.Jy:
        """Calculate the Degrees Per Flux Unit (DPFU) gain (also known as forward or telescope gain) of the antenna. A physical property of a telescope that describes how its 
        receivers respond to an increase in Janskys

        Returns
        -------
            (K/Jy)

        Reference
        ---------
            (1) https://casper.astro.berkeley.edu/astrobaki/index.php/Radiometer_Equation_Applied_to_Telescopes
        """
        if self._gain_dpfu is None:
            if self._is_set([self._effective_area]):
                # (1)
                self._gain_dpfu = (self._effective_area / (2 * joule_to_jansky(const.k_B)))
                self._gain_dpfu_eqn = """
G_f = \frac{A_e}{2k}      (m^2 / Jy)

where 
    $A_e$ is the antenna's effective area  ($m^2$)
    $k$ is Boltzmann's constant($J/s^{-1}/m^{-2}$)
"""
        return self._gain_dpfu

    @property
    def radio_gain(self) -> float:
        """Calculate the radio gain of the antenna.

        Returns
        -------
            (W)

        Reference
        ---------
            (1) https://uk.mathworks.com/help/phased/ref/aperture2gain.html (A. 12)
        """
        if not self._is_set([self._radio_gain]):
            if self._is_set([self._wavelength, self._effective_area]):
                # (3)
                self._radio_gain = ((self._effective_area  * 4.0 * pi) / pow(self._wavelength, 2))
                self._radio_gain_eqn = """
G_r = \frac{4 \pi A_e}{\lamba^2}      (W)

where
    $A_e$ is the antenna's effective area ($m^2$)
    $\lambda$ is the wavelength ($m$)
"""
        return self._radio_gain
    
    @property
    def effective_area(self) -> u.m * u.m:
        """Calculate the astronomy gain of the antenna.

        Returns
        -------
            (m^2)

        Reference:
        ---------
            (3) https://home.ifa.hawaii.edu/users/jpw/classes/radio/lectures/antenna_fundamentals.pdf
        """
        if self._effective_area is None:
            if self._is_set(
                [self.n_element_coherent, self.n_element_incoherent, self.centre_frequency]
            ):
                # (1)
                self._effective_area = (
                    self.n_element_coherent
                    * sqrt(self.n_element_incoherent)
                    * pow(frequency_to_wavelength(self.centre_frequency), 2)
                    / 2.0
                    / pi
                )
            if self._is_set([self.diameter]):
                # (2)
                self._effective_area = pi * pow(self.diameter / 2.0, 2)
        return self._effective_area

    def is_observable(self, pulsar: SkyCoord):
        return self._signal_to_noise(pulsar) > self.signal_to_noise

    @property
    def centre_frequency(self) -> u.Hz:
        if self._centre_frequency is None:
            if self._is_set([self._wavelength]):
                self._centre_frequency = wavelength_to_frequency(self._wavelength)
        return self._centre_frequency

    @property
    def wavelength(self) -> u.m:
        if self._wavelength is None:
            if self._is_set([self._centre_frequency]):
                self._wavelength = frequency_to_wavelength(self._centre_frequency)
        return self._wavelength

    @property
    def noise_power(self) -> u.W:
        """
        Reference: 
            https://www.researchgate.net/publication/280085456_Radio_Pulsar_Receiver_Systems_for_Space_Navigation
        """
        if self._noise_power is None:
            if self._is_set([self._system_temp, self._galaxy_noise_temp, self._solar_noise_temp]):
                self._noise_power = (
                        jansky_to_watt(const.k_B) * (
                        self._system_temp +
                        BACKGROUND_TEMP +
                        self._galaxy_noise_temp + 
                        self._solar_noise_temp
                    )
                    * self._bandwidth
                ).to(u.W)
        return self._noise_power
    
    def snr(self, pulsar: SkyCoord = None, integration_time: u.s = None) -> u.W:
        snr = None
        if self._snr is None:
            # snr = (Smean*jy) * Aeff * sqrt(np*deltaNu*tobs)/2.0/k/(tsys+tsky)*sqrt((psr[i].p0-w)/w);
            if self._is_set([pulsar, self._noise_power]):
                receive_power = self.receive_power(pulsar)
                snr = receive_power / self._noise_power
            elif self._is_set([pulsar, integration_time, self._n_polarisations, self._bandwidth, self._pulse_peak_amplitude, 
                self._system_temp]):
                """
                Reference:
                Buist, P.J., Engelen, S., Noroozi, A., Sundaramoorthy, P., Verhagen, S., Verhoeven, C., 
                "Principles and Potential of Pulsar Navigation," Proceedings of the 24th International 
                Technical Meeting of the Satellite Division of The Institute of Navigation (ION GNSS 2011), 
                Portland, OR, September 2011, pp. 3503-3515. 
                Equation  (15)
                """
                snr = (
                    sqrt(self._n_polarisations * integration_time * self._bandwidth) 
                    * (pulsar.ref_flux_density * self.gain_dpfu / self._system_temp)
                    * sqrt(pulsar.pulse_width(pulsar.period - pulsar.pulse_width))
                    / pulsar.period
                ) 
        else:
            snr = self._snr
        return snr

    def receive_power(self, pulsar: SkyCoord) -> u.W:
        if self._is_set([self._effective_area, self._polarisation_param, self._centre_frequency, self._bandwidth]):
            receive_power = (
                self._polarisation_param * self._effective_area 
                * jansky_to_watt(pulsar.flux_density(self._centre_frequency))
                * self._bandwidth
            ).to(u.W)
        return receive_power

    def toa_err(self, ref_antenna: "Antenna" = None, ref_integration_time: u.s = None, integration_time: u.s = None) -> u.ms:
        if self._toa_err is None:
            if self._is_set([self._bandwidth, self._snr, ref_antenna, ref_integration_time, integration_time, ref_antenna.bandwidth, ref_antenna.snr, ref_antenna.toa_err()]):
                toa_err = np.sqrt(
                    ref_antenna.toa_err()**2 * (ref_antenna.snr()**2 * ref_antenna.bandwidth * ref_integration_time) 
                    / (self.snr()**2 * self.bandwidth * integration_time) 
                )
            # Assume same integration time
            if self._is_set([self._bandwidth, self._snr, ref_antenna, ref_antenna.bandwidth, ref_antenna.snr(), ref_antenna.toa_err()]):
                toa_err = np.sqrt((
                    (ref_antenna.toa_err()**2 * ref_antenna.snr()**2 * ref_antenna.bandwidth)
                    / (self._snr**2 * self._bandwidth) 
                ))
        else:
            toa_err = self._toa_err
        return toa_err

    # @property
    # def timing_quality(self, pulsar: SkyCoord) -> u.dimensionless_unscaled:
    #     q_timing = sum([const.k_B / (self.system_temp * )])
        
    
    @property
    def galaxy_noise_temp(self) -> u.K:
        """
        background noise temperature from the galaxy

        Reference: 
            https://www.researchgate.net/publication/280085456_Radio_Pulsar_Receiver_Systems_for_Space_Navigation
        """
        if self._galaxy_noise_temp is None:
            self._galaxy_noise_temp = 6 * (self.centre_frequency.to_value(u.GHz))**(-2.2) * u.K
        return self._galaxy_noise_temp

    @property
    def solar_noise_temp(self) -> u.K:
        """
        solar system noise temperature

        Reference: 
            https://www.researchgate.net/publication/280085456_Radio_Pulsar_Receiver_Systems_for_Space_Navigation

        TODO fix units
        """
        if self._solar_noise_temp is None:
            if self._is_set([self._centre_frequency, self._effective_area, self._side_lobe_attenuation, self._solar_distance]):
                self._solar_noise_temp = (
                    (72 * self._centre_frequency.to(u.GHz) + 0.058 * u.GHz) 
                    * self._effective_area 
                    * 10**(self._side_lobe_attenuation / 10)
                    * self._solar_distance.to_value(u.AU)**(-2)
                )
        return self._solar_noise_temp 

    @property
    def polsarisation_param(self) -> int:
        if self._polsarisation_param is None:
            if self._is_set([self._n_polarisations,]) and 0 <= self._n_polarisations <= 2:
                self._polarisation_param = self._n_polarisations / 2
        return self._polsarisation_param
    
    @property
    def n_element_coherent(self) -> int:
        return self._n_element_coherent
    
    @property
    def n_element_incoherent(self) -> int:
        return self._n_element_incoherent
    
    @property
    def n_polarisations(self) -> int:
        return self._n_polarisations

    @property
    def system_temp(self) -> u.K:
        return self._system_temp
    
    @property
    def side_lobe_attenuation(self) -> float:
        return self._side_lobe_attenuation

    @property
    def solar_distance(self) -> u.AU:
        return self._solar_distance

    @property
    def bandwidth(self) -> u.MHz:
        return self._bandwidth
    
    @property
    def diameter(self) -> u.m:
        return self._diameter
    
    @property
    def elevation(self) -> u.m:
        return self._elevation
    
    @property
    def location(self) -> u.m:
        return self._location
    
    
    @property
    def radio_gain_equation(self) -> str:
        return self._radio_gain_eqn
    
    @property
    def gain_dpfu_equation(self) -> str:
        return self._gain_dpfu_eqn
    
    @property
    def centre_frequency(self) -> u.MHz:
        return self._centre_frequency

    def _is_set(self, vals: list[Any]) -> bool:
        return None not in vals
    
    def _propagate_calculations(self) -> float:
        for _ in range(3):
            for var in type(self).__dict__:
                getattr(self, var)

