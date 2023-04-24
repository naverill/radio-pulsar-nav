from math import pi

from astropy.coordinates import EarthLocation

from radio_pulsar_nav.antenna import Antenna

kens = Antenna(
    name="Marsfield telescope",
    signal_to_noise=1,
    temp=60.0,  # K
    effective_area=pi * 2**2,
    bandwidth=40e6,  # (Hz)
    centre_freq=1400,  # (MHz)
    pos=EarthLocation.from_geodetic(lat=-33.77290643916046, lon=151.0976937264337),
)

parkes = Antenna(
    name="Parkes telescope",
    signal_to_noise=10,  #
    temp=25.0,  # K
    gain=0.6,  # (K/Jy)
    bandwidth=340e6,  # (Hz)
    centre_freq=1374,  # (MHz)
    pos=EarthLocation.from_geodetic(lat=-32.99327814611731, lon=148.26503125664433),
)

fast = Antenna(
    name="FAST Telescope",
    signal_to_noise=9.0,  #
    temp=25.0,  # (K)
    gain=16.0,  # (K/Jy)
    bandwidth=512e6,  # (Hz)
    centre_freq=1350,  # (MHz)
    pos=EarthLocation.from_geodetic(lat=25.654006939684034, lon=106.85784898726897),
)
