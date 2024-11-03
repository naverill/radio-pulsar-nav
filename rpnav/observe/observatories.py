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
    name="Parkes",
    system_temp=21.0 * u.K,  # K
    gain_dpfu=(1.0 / 1.8) * u.K / u.Jy,  # (K/Jy)
    time=Time(60002.3, format="mjd"),
    bandwidth=(128 * u.MHz).to(u.Hz),  # (Hz)
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

NAVIGATE = Antenna(
        name="navigate",
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

WOODCHESTER = Antenna(
        name="woodchester",
        system_temp=30.0 * u.K,  # K
        bandwidth=(8 * u.MHz).to(u.Hz),  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(820 * u.MHz).to(u.Hz),  # (MHz)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),
        diameter=2.4 * u.m
)

"""
Reference:
    The second data release from the European Pulsar Timing Array
    LEAP: the Large European Array for Pulsars

The typical integration length per source was 30 min
"""
EFFELSBERG = Antenna(
        name="effelsberg",
        bandwidth=(500 * u.MHz).to(u.Hz),  # (Hz) 
        system_temp=24.0 * u.K,  # K
        gain_dpfu=1.5 * u.K / u.Jy,  # (K/Jy)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(1400 * u.MHz).to(u.Hz),  # (MHz)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),         # TODO fix
        diameter=100 * u.m
)

"""
Reference:
    The second data release from the European Pulsar Timing Array
    LEAP: the Large European Array for Pulsars

The typical integration time per source was varied depend-
ing on the pulsar and epoch, with the median observation times
per pulsar ranging between ∼10 min to 55 min
"""
LOVELL = Antenna(
        name="lovell",
        bandwidth=(384 * u.MHz).to(u.Hz),  # (Hz)
        time=Time(60002.3, format="mjd"),
        system_temp=25.0 * u.K,  # K
        gain_dpfu=1 * u.K / u.Jy,  # (K/Jy)
        centre_frequency=(1.4 * u.MHz).to(u.Hz),  # (MHz)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),  # TODO fix
        diameter=76 * u.m
)

"""
Reference:
    The second data release from the European Pulsar Timing Array
    LEAP: the Large European Array for Pulsars

NUPPI observations have durations ranging from
20 min to 80 min, and cover a frequency bandwidth of 512 MHz
channelised into 128 channels that are coherently dedispersed in
real-time. Observations with the L-band receiver were made at
a central frequency of 1484 MHz. Those with the S -band re-
ceiver were generally centred on 2539 MHz and occasionally
on 1854 and 2154 MHz
"""
NRT = Antenna(
        name="NRT",
        bandwidth=(512 * u.MHz).to(u.Hz),  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(1484 * u.MHz).to(u.Hz),  # (MHz)
        system_temp=35 * u.K,  # K
        gain_dpfu=1.4 * u.K / u.Jy,  # (K/Jy)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),  # TODO fix
        diameter=94 * u.m
)

"""
Reference:
    The second data release from the European Pulsar Timing Array
    LEAP: the Large European Array for Pulsars
"""
SRT = Antenna(
        name="SRT",
        bandwidth=(500 * u.MHz).to(u.Hz),  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(1484 * u.MHz).to(u.Hz),  # (MHz)
        system_temp=20 * u.K,  # K
        gain_dpfu=0.63 * u.K / u.Jy,  # (K/Jy)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),  # TODO fix
        diameter=64 * u.m
)

"""
Reference:
    The second data release from the European Pulsar Timing Array
    LEAP: the Large European Array for Pulsars

The duration of pulsar observations at the WSRT varied with
pulsar and epoch and lasted between 15 min and 45 min
"""
WSRT = Antenna(
        name="WSRT",
        bandwidth=(160 * u.MHz).to(u.Hz),  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(1380 * u.MHz).to(u.Hz),  # (MHz) 
        system_temp=27 * u.K,  # K
        gain_dpfu=1.2 * u.K / u.Jy,  # (K/Jy)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),  # TODO fix
        diameter=94 * u.m
)

"""
Reference:
    The second data release from the European Pulsar Timing Array

The typical exposure time per pair of pulsar and phase
calibrator is 45–60 and 2–3 min

The LEAP S/N should then be similar to the sum of the S/Ns of the individual telescopes
The peak S/Ns of Effelsberg, Jodrell Bank, Nanc¸ay, WSRT and LEAP are
97, 51, 42, 30 and 220,
"""
LEAP = Antenna(
        name="LEAP",
        # gain_dpfu=(1.0 / 1.8) * u.K / u.Jy,  # (K/Jy)
        bandwidth=(128 * u.MHz).to(u.Hz),  # (Hz)
        time=Time(60002.3, format="mjd"),
        centre_frequency=(1400 * u.MHz).to(u.Hz),  # (MHz)
        location=EarthLocation.from_geocentric(
            -4554231.5 * u.m,
            2816759.1 * u.m, 
            -3454036.3 * u.m
        ),  # TODO fix
        diameter=194 * u.m
)

"""
Reference antenna defined in Feasibility Study (Sala, 2004)
"""
SALA = Antenna(
    name="SALA",
    centre_frequency=1 * u.GHz,
    bandwidth=200 * u.MHz,
    effective_area=10 * u.m * u.m,
    location=PARKES.location,
    time=PARKES.time,
)