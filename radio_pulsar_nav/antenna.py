"""
Antenna observer object
"""
from math import pi

import astropy.units as u
import numpy as np
from astroplan import Observer
from astropy.coordinates import EarthLocation

from radio_pulsar_nav.constants import BOLTZMANN_CONSTANT, SPEED_OF_LIGHT


class Antenna(Observer):
    """
    Antenna observer object
    """

    def __init__(
        self,
        name: str,
        signal_to_noise: float,
        temp: float,
        bandwidth: float,
        centre_freq: float,
        pos: EarthLocation,
        effective_area: float = None,
        gain: float = None,
    ):
        """Instantiate antenna.

        Parameters
        ----------
            signal_to_noise (float): Signal to noise ratio
            temp (float, K): System temperature
            gain (float, K/Jy): Antenna gain
            bandwidth (float, Hz): Antenna bandwidth
            centre_freq (float, MHz):
            half_beamwidth (float, degrees):
        """
        self.signal_to_noise = signal_to_noise
        self.temp = temp
        self.bandwidth = bandwidth
        self.centre_freq = centre_freq
        self._gain = gain
        self.effective_area = effective_area
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
        """
        return (
            (correction * self.signal_to_noise * self.temp)
            / (self.gain * np.sqrt(num_polarisations * self.bandwidth * integration_time))
            * np.sqrt(pulse_width_to_period / (1 - pulse_width_to_period))
        ) * 1000

    @property
    def gain(self):
        """Calculate the astronomy gain of the antenna.

        Returns
        -------
            (K/Jy)
        """
        if self._gain is None:
            self._gain = 2 * BOLTZMANN_CONSTANT / self.effective_area
        return self._gain
