import numpy as np

class Antenna:
    def __init__(
        self,
        signal_to_noise: float,
        temp: float,
        gain: float,
        bandwidth: float,
    ):
        self.signal_to_noise = signal_to_noise
        self.temp = temp
        self.gain = gain
        self.bandwidth = bandwidth

    def min_observable_flux_density(
        self,
        integration_time: float,
        pulse_width: float,
        pulse_period: float,
        correction: float = 1,
        num_polarisations: float = 2,
    ):
        """
        Parameters
        ----------
        """
        return (
            (correction * signal_to_noise * temp)
            / (
                gain * np.sqrt(num_polarisations * bandwidth * integration_time)
            )
            * np.sqrt(
                pulse_width / (period - pulse_width)
            )
        )
