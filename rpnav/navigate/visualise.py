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


def plot_scattermap(long: list[float], lat: list[float], err: list[float]):
    fig = go.Figure(
        go.Scattermap(
            # mode = "markers+lines",
            mode = "markers",
            lon = long,
            lat = lat,
            # marker = {'size': 1},
            line={'width':1},
            opacity=1,
        )
    )
    fig.update_layout(
        margin ={'l':0,'t':0,'b':0,'r':0},
        width=1800,
        height=1800,
        map = {
            # 'style': "white-bg",
            'style': "open-street-map",
            'center': {'lon': 148.26356, 'lat': -32.99841},
            'zoom': 2})
    return fig

