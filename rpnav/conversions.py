from rpnav.constants import SPEED_OF_LIGHT


def frequency_to_wavelength(freq: float) -> float:
    """
    Convert frequency to wavelength 

    Args:
        freq (float, Hz):

    Returns:
        Wavelength (m)
    """
    return SPEED_OF_LIGHT / freq


def wavelength_to_frequency(wav: float) -> float:
    """
    Convert wavelength to frequency 
    
    Args:
        wavelength (m)

    Returns:
        Frequency (Hz)
    """
    return SPEED_OF_LIGHT / wav
