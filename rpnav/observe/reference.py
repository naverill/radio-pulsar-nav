import numpy as np
import astropy.units as u

from rpnav.observe.antenna import Antenna
from rpnav.observe.observation import Observation
from rpnav.observe.pulsar import Pulsar
from rpnav.observe.observatories import PARKES
from rpnav.observe.observatories import SALA
from rpnav.observe.observatories import NRT
from rpnav.observe.observatories import LEAP
from rpnav.observe.observatories import WSRT
from astropy import constants as const


"""
Define set of known reference observations and use as a reference with 

References:
    The second data release from the European Pulsar Timing Array
"""
OBSERVATIONS = []
pulsars = Pulsar.load_catalogue()
for p in pulsars:
    """
    Reference:
        An ultra-wide bandwidth (704 to 4 032 MHz) receiver for the Parkes radio telescope (Hobbs, 2020)
    """
    if p.jname == "J1909-3744":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(145 * u.ns),
                integration_time= 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        """
        Reference:
            Timing analysis for 20 millisecond pulsars in the Parkes Pulsar Timing Array (Reardon, 2016)
        """
    elif p.jname == "J0711-6830":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(2.0 * u.us).to(u.ns),
                integration_time= 556 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
    elif p.jname == "J1024-0719":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(10.4 * u.us).to(u.ns),
                integration_time= 493 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (220 * u.dB(u.dimensionless_unscaled)),
                toa_err=(2.38 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J1730-2304":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(1.9 * u.us).to(u.ns),
                integration_time= 390 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(1.57 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J1744-1134":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(0.5 * u.us).to(u.ns),
                integration_time= 534 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (220 * u.dB(u.dimensionless_unscaled)),
                toa_err=(0.90 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname in ["J1824-2452A", "B1821-24A"]:
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(5.5 * u.us).to(u.ns),
                integration_time= 313 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
    elif p.name in ["J1939+2134", "B1937+21"]:
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(5.8 * u.us).to(u.ns),
                integration_time= 397 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                SALA,
                p,
                snr= (-55.6 * u.dB(u.dimensionless_unscaled)),
                toa_err=(3162 * u.m / const.c).to(u.ns),
                integration_time=(36.36 * u.min).to(u.s),
                name=f"{p.jname}_{SALA.name}",
            )
        )
    elif p.jname == "J2124-3358":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(2.9 * u.us).to(u.ns),
                integration_time= 652 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(3.7 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J0613-0200":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(1.0 * u.us).to(u.ns),
                integration_time= 639 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (220 * u.dB(u.dimensionless_unscaled)),
                toa_err=(1.43 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J1045-4509":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(11.9 * u.us).to(u.ns),
                integration_time= 646 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
    elif p.jname == "J1603-7202":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(3.1 * u.us).to(u.ns),
                integration_time= 244 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
    elif p.jname in ["J1857+0943", "B1855+09"]:
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(1.1 * u.us).to(u.ns),
                integration_time= 291 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (220 * u.dB(u.dimensionless_unscaled)),
                toa_err=(1.70 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J2129-5721":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(1.4 * u.us).to(u.ns),
                integration_time=448 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
    elif p.jname == "J2145-0750":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(1.7 * u.us).to(u.ns),
                integration_time=972 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
    elif p.jname == "J1022+1001":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(1.8 * u.us).to(u.ns),
                integration_time=615 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (220 * u.dB(u.dimensionless_unscaled)),
                toa_err=(2.19 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J1600+3053":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(0.8 * u.us).to(u.ns),
                integration_time=715 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (220 * u.dB(u.dimensionless_unscaled)),
                toa_err=(0.48 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J1643-1224":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(2.8 * u.us).to(u.ns),
                integration_time=488 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
    elif p.jname == "J0437-4715":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(0.3 * u.us).to(u.ns),
                integration_time=5065 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
    elif p.jname == "J1713+0747":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(0.4 * u.us).to(u.ns),
                integration_time=622 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (220 * u.dB(u.dimensionless_unscaled)),
                toa_err=(0.32 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J1909-3744":
        OBSERVATIONS.append(
            Observation(
                PARKES,
                p,
                snr=146.09 * u.dimensionless_unscaled,
                toa_err=(0.2 * u.us).to(u.ns),
                integration_time=1368 * 1800 * u.s,
                name=f"{p.jname}_{PARKES.name}",
            )
        )
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(0.29 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    # elif p.name in ["J0835-4510", "J1017-7156", "J1732-5049", "J2241-5236", "J1125-6014", "J1446-4701", "J1545-4550", "J1713+0747", "J1832-0836"]:
    #     OBSERVATIONS.append(
    #         Observation(
    #             PARKES,
    #             p,
    #             snr=146.09 * u.dimensionless_unscaled,
    #             toa_err=(145 * u.ns).to(u.ns), # TODO fix value
    #             integration_time= 1800 * u.s,
    #             name=f"{p.jname}_{PARKES.name}",
    #         )
    #     )
        """
        Reference:
            Feasibility Study for a Spacecraft Navigation System relying on Pulsar Timing Information Final Report (Sala, 2004)
        """
    elif p.name in ["B0736-40", "J0738-4042"]: 
        OBSERVATIONS.append(
            Observation(
                SALA,
                p,
                snr= (-50.2 * u.dB(u.dimensionless_unscaled)),
                toa_err=(9246467 * u.m / const.c).to(u.ns),
                integration_time=(3.07 * u.min).to(u.s),
                name=f"{p.jname}_{SALA.name}",
            )
        )
    elif p.name in ["B1451-68", "J1456-6843"]:
        OBSERVATIONS.append(
            Observation(
                SALA,
                p,
                snr= (-50.0 * u.dB(u.dimensionless_unscaled)),
                toa_err=(452529 * u.m / const.c).to(u.ns),
                integration_time=(2.35 * u.min).to(u.s),
                name=f"{p.jname}_{SALA.name}",
            )
        )
    elif p.name == "B0950-08": 
        OBSERVATIONS.append(
            Observation(
                SALA,
                p,
                snr= (-49.8 * u.dB(u.dimensionless_unscaled)),
                toa_err=(316174 * u.m / const.c).to(u.ns),
                integration_time=(1.04 * u.min).to(u.s),
                name=f"{p.jname}_{SALA.name}",
            )
        )
    elif p.name in ["B0329+54", "J0332+5434"]: 
        OBSERVATIONS.append(
            Observation(
                observer=SALA,
                pulsar=p,
                snr= (-45.2 * u.dB(u.dimensionless_unscaled)),
                toa_err=(290556 * u.m / const.c).to(u.ns),
                integration_time=(0.12 * u.min).to(u.s),
                name=f"{p.jname}_{SALA.name}",
            )
        )
    elif p.jname == "J0030+0450": 
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (125 * u.dB(u.dimensionless_unscaled)),
                toa_err=(3.40 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J0751-1807": 
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (125 * u.dB(u.dimensionless_unscaled)),
                toa_err=(2.22 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J0900-3144": 
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(2.95 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J1012+5307": 
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (125 * u.dB(u.dimensionless_unscaled)),
                toa_err=(1.76 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J1455-3330": 
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(7.23 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J1640+2224": 
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (125 * u.dB(u.dimensionless_unscaled)),
                toa_err=(3.57 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J1738+0333": 
        OBSERVATIONS.append(
            Observation(
                WSRT,
                p,
                snr= (30 * u.dB(u.dimensionless_unscaled)),
                toa_err=(4.31 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{WSRT.name}",
            )
        )
    elif p.jname == "J1751+2857": 
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(3.17 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J1801-1417": 
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(4.09 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J1804-2717": 
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(5.94 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J143-1113": 
        OBSERVATIONS.append(
            Observation(
                WSRT,
                p,
                snr= (30 * u.dB(u.dimensionless_unscaled)),
                toa_err=(1.37 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{WSRT.name}",
            )
        )
    elif p.jname == "J1910+1256": 
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(2.65 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J1911+1347": 
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr= (42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(1.22 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
    elif p.jname == "J1918+0642": 
        OBSERVATIONS.append(
            Observation(
                LEAP,
                p,
                snr= (220 * u.dB(u.dimensionless_unscaled)),
                toa_err=(2.04 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{LEAP.name}",
            )
        )
    elif p.jname == "J2322+2057": 
        OBSERVATIONS.append(
            Observation(
                NRT,
                p,
                snr=(42 * u.dB(u.dimensionless_unscaled)),
                toa_err=(10.9 * u.us).to(u.ns),
                integration_time=(60 * u.min).to(u.s),
                name=f"{p.jname}_{NRT.name}",
            )
        )
