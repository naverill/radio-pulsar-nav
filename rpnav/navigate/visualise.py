from math import log2

import astropy.units as u
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from astropy.coordinates import ICRS, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time

from rpnav.observe.antenna import Antenna
from rpnav.pulsar import Pulsar

def plot_heatmap(rmseMap: list[list[float]]):
    fig = go.Figure(data=go.Heatmap(z=rmseMap),
        #     color__scale='RdBu_r', 
    )
    return fig

