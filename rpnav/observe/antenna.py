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

from rpnav.constants import BOLTZMANN_CONSTANT
from rpnav.conversions import frequency_to_wavelength, wavelength_to_frequency, joules_to_jansky
from rpnav.observe.observer import Observer


class Antenna(Observer):
    """
    Antenna observer object
    """

    def __init__(
        self,
        name: str,
        system_temp:  u.K,
        bandwidth: u.Hz,
        location: EarthLocation,
        snr: u.dimensionless_unscaled = None,
        aliases: list[str] = [],
        time: Time = None,
        timescale: str = "TCB",
        itoa_code: str = None,
        centre_frequency: u.MHz = None,
        wavelength: u.m = None,
        effective_area: u.m * u.m = None,
        diameter: u.m = None,
        n_polarisations: int = 2,
        n_element_coherent: int = None,
        n_element_incoherent: int = None,
        gain_dpfu: u.K / u.Jy = None,
        astronomy_gain: u.m * u.m / u.Jy = None,
        elevation: u = u.m,
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
        self._effective_area = effective_area
        self._n_element_coherent = n_element_coherent
        self._n_element_incoherent = n_element_incoherent
        self._diameter  : u.m = diameter
        self._radio_gain = None
        self._n_polarisations: int = n_polarisations
        self._elevation: u.m = elevation
        self._location : EarthLocation = location
        self._gain_dpfu_eqn = ""
        self._radio_gain_eqn = ""
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
            (correction * self.system_temp)
            / (np.sqrt(self.n_polarisations * self.bandwidth * integration_time) * self.gain_dpfu)
            * np.sqrt(pulse_width_to_period / (1 - pulse_width_to_period))
            * self.snr
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
                self._gain_dpfu = (self._effective_area / (2 * joules_to_jansky(BOLTZMANN_CONSTANT)))
                self._gain_dpfu_eqn = """
G_f = \frac{A_e}{2k}      (m^2 / Jy)

where 
    $A_e$ is the antenna's effective area  ($m^2$)
    $k$ is Boltzmann's constant($J/s^{-1}/m^{-2}$)
"""
        return self._gain_dpfu
    
    @property
    def gain_dpfu_equation(self) -> str:
        return self._gain_dpfu_eqn

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
    def radio_gain_equation(self) -> str:
        return self._radio_gain_eqn
    
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

    def pulsar_snr(self, pulsar: SkyCoord):
        # snr = (Smean*jy) * Aeff * sqrt(np*deltaNu*tobs)/2.0/k/(tsys+tsky)*sqrt((psr[i].p0-w)/w);
        raise NotImplementedError

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
    def centre_frequency(self) -> EarthLocation:
        return self._centre_frequency
    
    @property
    def system_temp(self) -> u.K:
        return self._system_temp
    
    @property
    def snr(self) -> u.W:
        return self._snr
    
    @property
    def n_element_coherent(self) -> int:
        return self._n_element_coherent
    
    @property
    def n_element_incoherent(self) -> int:
        return self._n_element_incoherent
    
    @property
    def n_polarisations(self) -> int:
        return self._n_polarisations

    def _is_set(self, vals: list[Any]) -> bool:
        return None not in vals
    
    def _propagate_calculations(self) -> float:
        for _ in range(3):
            for var in type(self).__dict__:
                getattr(self, var)

