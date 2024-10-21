import astropy.units as u
from astropy.constants import Constant

BOLTZMANN_CONSTANT = Constant(
    abbrev="k",
    name="Boltzmann Constant",
    value=1.380649e-23,
    uncertainty=0,
    reference="",
    unit= u.J / u.K # J/s-1/m-2
)



"""
1 Jansky = 10-26 Joules s-1 m-2 = 10-26 Watts m-2
"""