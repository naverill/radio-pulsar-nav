from math import pi

from astropy.coordinates import EarthLocation

from radio_pulsar_nav.antenna import Antenna

kens = Antenna(
    name="Marsfield telescope",
    signal_to_noise=1,
    system_temp=30.0,  # K
    bandwidth=40e6,  # (Hz)
    centre_frequency=1400,  # (MHz)
    pos=EarthLocation.from_geodetic(lat=-33.77290643916046, lon=151.0976937264337),
    diameter=2,
    # n_element_coherent = 1,
    # n_element_incoherent = 1,
)
parkes = Antenna(
    name="Parkes telescope",
    signal_to_noise=10,  #
    system_temp=25.0,  # K
    astronomy_gain=0.6,  # (K/Jy)
    bandwidth=340e6,  # (Hz)
    centre_frequency=1374,  # (MHz)
    pos=EarthLocation.from_geodetic(lat=-32.99327814611731, lon=148.26503125664433),
)

fast = Antenna(
    name="FAST Telescope",
    signal_to_noise=9.0,  #
    system_temp=25.0,  # (K)
    astronomy_gain=16.0,  # (K/Jy)
    bandwidth=512e6,  # (Hz)
    centre_frequency=1350,  # (MHz)
    pos=EarthLocation.from_geodetic(lat=25.654006939684034, lon=106.85784898726897),
)
