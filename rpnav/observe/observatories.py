from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u

from rpnav.observe.antenna import Antenna

"""
@brief Parkes ultra-wide-bandwidth low-frequency receiver (UWL) telescope config

@details 
    References: 
        https://www.parkes.atnf.csiro.au/observing/documentation/users_guide/pkug.html
        https://www.cambridge.org/core/journals/publications-of-the-astronomical-society-of-australia/article/an-ultrawide-bandwidth-704-to-4-032-mhz-receiver-for-the-parkes-radio-telescope/54E3119D8DFAC9C8D70635CE3C7E8099
"""
PARKES = Antenna(
        name="Parkes telescope",
        snr=10 * u.dimensionless_unscaled,  #
        system_temp=21.0 * u.K,  # K
        gain_dpfu=(1.0 / 1.8) * u.K / u.Jy,  # (K/Jy)
        time=Time(60002.3, format="mjd"),
        bandwidth=((4032 - 704) * u.MHz).to(u.Hz),  # (Hz)
        centre_frequency=(1400 * u.MHz).to(u.Hz),  # (MHz)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),
        diameter=64 * u.m
)

FAST = Antenna(
        name="FAST Telescope",
        snr=9.0 * u.dimensionless_unscaled,  #
        system_temp=25.0 * u.K,  # (K)
        gain_dpfu=16.0 * u.K / u.Jy,  # (K/Jy)
        bandwidth=(3300 - 500) * u.MHz,  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=1350 * u.MHz,  # (MHz)
        location=EarthLocation.from_geodetic(
            lat=25.654006939684034 * u.deg, lon=106.85784898726897 * u.deg
        ),
        diameter=500 * u.m
    )


MSFD = Antenna(
        name="msfd",
        snr=1 * u.dimensionless_unscaled,
        system_temp=30.0 * u.K,  # K
        bandwidth=40e6 * u.Hz,  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(1400 * u.MHz).to(u.Hz),  # (MHz)
        location=EarthLocation.from_geocentric(
            x=-4646251.90 * u.m,
            y=2565112.87 * u.m,
            z=-3525535.83 * u.m,
        ),
        diameter=2 * u.m
    )

WOODCHESTER_WEAK = Antenna(
        name="woodchester_weak",
        snr= 10.2797461928934 * u.dimensionless_unscaled,  #
        # gain_dpfu=(1.0 / 1.8) * u.K / u.Jy,  # (K/Jy)
        system_temp=30.0 * u.K,  # K
        bandwidth=40e6 * u.Hz,  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(840 * u.MHz).to(u.Hz),  # (MHz)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),
        diameter=2.4 * u.m
)

WOODCHESTER_STRONG = Antenna(
        name="woodchester_strong",
        snr=15.17340361445783 * u.dimensionless_unscaled,  #
        # gain_dpfu=(1.0 / 1.8) * u.K / u.Jy,  # (K/Jy)
        system_temp=30.0 * u.K,  # K
        bandwidth=40e6 * u.Hz,  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(840 * u.MHz).to(u.Hz),  # (MHz)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),
        diameter=2.4 * u.m
)