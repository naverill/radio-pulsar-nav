from rpnav.constants import SPEED_OF_LIGHT


def frequency_to_wavelength(freq: float):
    """
    Parameters
    ----------
        freq (float, Hz):
    """
    return SPEED_OF_LIGHT / freq


def wavelength_to_frequency(wav: float):
    """
    Parameters
    ----------
    wav (float, m)
    """
    return SPEED_OF_LIGHT / wav
