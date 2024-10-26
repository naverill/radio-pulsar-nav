import astropy.units as u
from astropy.constants import Constant

BACKGROUND_TEMP = Constant(
    abbrev="k",
    name="Background Temperature",
    value=2.7,
    uncertainty=0,
    reference="",
    unit= u.K # J/s-1/m-2
)


"""
1 Jansky = 10-26 Joules s-1 m-2 = 10-26 Watts m-2
"""