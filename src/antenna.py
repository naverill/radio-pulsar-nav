import numpy as np

class Antenna:
    def __init__(
        self,
        name: str,
        signal_to_noise: float,
        temp: float,
        gain: float,
        bandwidth: float,
        centre_freq: float,
        half_beamwidth: float,
        lat: float,
        long: float,
        min_abs_galactic_lat: float,
        max_abs_galactic_lat: float,
    ):
        """
        Parameters
        ----------
            signal_to_noise (float): Signal to noise ratio
            temp (float, K): System temperature
            gain (float): Antenna gain
            bandwidth (float, Hz): Antenna bandwidth
            centre_freq (float, MHz): 
            half_beamwidth (float, degrees):
        """
        self.name = name
        self.signal_to_noise = signal_to_noise
        self.temp = temp
        self.gain = gain
        self.bandwidth = bandwidth
        self.centre_freq = centre_freq
        self.half_beamwidth = half_beamwidth
        self.lat = lat
        self.long = long
        self.min_abs_galactic_lat = min_abs_galactic_lat
        self.max_abs_galactic_lat = max_abs_galactic_lat

    def min_observable_flux_density(
        self,
        integration_time: float,
        pulse_width_to_period: float = 0.1,
        correction: float = 1,
        num_polarisations: int = 2,
    ):
        """
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
            / (
                self.gain * np.sqrt(num_polarisations * self.bandwidth * integration_time)
            )
            * np.sqrt(
                pulse_width_to_period / (1 - pulse_width_to_period)
            )
        ) * 1000
