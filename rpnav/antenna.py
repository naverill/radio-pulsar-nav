"""
Antenna observer object
"""
from math import pi, sqrt
from typing import Any

import astropy.units as u
import numpy as np
from astroplan import Observer
from astropy.coordinates import EarthLocation, SkyCoord

from rpnav.constants import BOLTZMANN_CONSTANT, SPEED_OF_LIGHT
from rpnav.conversions import frequency_to_wavelength, wavelength_to_frequency


class Antenna(Observer):
    """
    Antenna observer object
    """

    def __init__(
        self,
        name: str,
        signal_to_noise: float,
        system_temp: float,
        bandwidth: float,
        pos: EarthLocation,
        centre_frequency: float = None,
        wavelength: float = None,
        effective_area: float = None,
        diameter: float = None,
        n_element_coherent: int = None,
        n_element_incoherent: int = None,
        astronomy_gain: float = None,
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
        self.signal_to_noise = signal_to_noise
        self.system_temp = system_temp
        self.bandwidth = bandwidth
        self._centre_frequency = centre_frequency
        self._wavelength = wavelength
        self._astronomy_gain = astronomy_gain
        self._effective_area = effective_area
        self.n_element_coherent = n_element_coherent
        self.n_element_incoherent = n_element_incoherent
        self.diameter = diameter
        self._radio_gain = None
        self._propagate_calculations()
        super().__init__(name=name, location=pos, elevation=0 * u.m)

    def min_observable_flux_density(
        self,
        integration_time: float,
        pulse_width_to_period: float = 0.1,
        correction: float = 1,
        num_polarisations: int = 2,
    ):
        """Calculate minimum flux density that can be measured by the antenna
        given a fixed integration time.

        Parameters
        ----------
            integration_time (float, s):
            pulse_width_to_period (float):
            correction (float):
            num_polarisations (int):

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
            / (np.sqrt(num_polarisations * self.bandwidth * integration_time) * self.astronomy_gain)
            * np.sqrt(pulse_width_to_period / (1 - pulse_width_to_period))
            * self.signal_to_noise
            * 1000
        )

    @property
    def astronomy_gain(self) -> float:
        """Calculate the astronomy gain of the antenna.

        Returns
        -------
            (K/Jy)

        Reference
        ---------
            (1) https://link.springer.com/content/pdf/bbm:978-3-642-19627-0/1.pdf (A. 12)
        """
        if self._astronomy_gain is None:
            if self._is_set([self._effective_area]):
                # (1)
                self._astronomy_gain = self._effective_area / (2 * BOLTZMANN_CONSTANT) * 1e-26
        return self._astronomy_gain

    @property
    def radio_gain(self) -> float:
        """Calculate the astronomy gain of the antenna.

        Returns
        -------
            (K/Jy)

        Reference
        ---------
            (1) https://link.springer.com/content/pdf/bbm:978-3-642-19627-0/1.pdf (A. 12)
        """
        if not self._is_set([self._radio_gain]):
            if self._is_set([self._wavelength, self._effective_area]):
                # (3)
                self._radio_gain = self._effective_area * pow(self._wavelength, 2) / 4.0 / pi
        return self._radio_gain

    @property
    def effective_area(self) -> float:
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
                    * pow(SPEED_OF_LIGHT / self.centre_frequency * 1e6, 2)
                    / 2.0
                    / pi
                )
            if self._is_set([self.diameter]):
                # (2)
                self._effective_area = pi * pow(self.diameter / 2.0, 2)
        return self._effective_area

    def _is_set(self, vals: list[Any]) -> bool:
        return None not in vals

    def _signal_to_noise(self, pulsar: SkyCoord):
        # snr = (Smean*jy) * Aeff * sqrt(np*deltaNu*tobs)/2.0/k/(tsys+tsky)*sqrt((psr[i].p0-w)/w);
        raise NotImplementedError

    def is_observable(self, pulsar: SkyCoord):
        return self._signal_to_noise(pulsar) > self.signal_to_noise

    @property
    def centre_frequency(self) -> float:
        if self._centre_frequency is None:
            if self._is_set([self._wavelength]):
                self._centre_frequency = wavelength_to_frequency(self._wavelength)
        return self._centre_frequency

    @property
    def wavelength(self) -> float:
        if self._wavelength is None:
            if self._is_set([self._centre_frequency]):
                self._wavelength = frequency_to_wavelength(self._centre_frequency)
        return self._wavelength

    def _propagate_calculations(self) -> float:
        for _ in range(3):
            for var in type(self).__dict__:
                getattr(self, var)
