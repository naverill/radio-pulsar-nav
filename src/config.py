from antenna import Antenna

parkes = Antenna(
    name="Parkes telescope",
    signal_to_noise = 10,   # 
    temp = 25.,             # K
    gain = 0.6,             # 
    bandwidth = 340e6,      # (Hz)
    centre_freq = 1374,     # (MHz)
    half_beamwidth = 7. * 0.0166667,    # (degrees)
    long = -32.99327814611731,
    lat = 148.26503125664433,
    min_abs_galactic_lat = 20,
    max_abs_galactic_lat = 45,
)

fast = Antenna(
    name="FAST Telescope",
    signal_to_noise = 9.0,  #
    temp = 25.,             # (K)
    gain = 16.,             #
    bandwidth = 512e6,      # (Hz)
    centre_freq = 1350,     # (MHz)
    half_beamwidth = 1. * 0.0166667,    # (degrees)
    lat=106.85784898726897,
    long=25.654006939684034,
    min_abs_galactic_lat = 0,
    max_abs_galactic_lat = 90,
)
