import astropy.units as u
from astropy import constants as const

def frequency_to_wavelength(value: u.Hz) -> u.m:
    """
    Convert frequency to wavelength 

    Args:
        freq (float, Hz):

    Returns:
        Wavelength (m)
    """
    return (const.c / value.to(1 / u.s)).to(u.m)


def wavelength_to_frequency(value: u.m) -> u.Hz:
    """
    Convert wavelength to frequency 
    
    Args:
        wavelength (m)

    Returns:
        Frequency (Hz)
    """
    return (const.c / value).to(u.Hz)

def joules_to_jansky(value: u.J) -> u.Jy * u.m * u.m:
    """
    Convert Joules to Jansky

    Reference:
        https://lweb.cfa.harvard.edu/~dfabricant/huchra/ay145/units.html
    """
    return (value * 1e26 * u.Jy * u.m * u.m # * u.s
            / u.J)