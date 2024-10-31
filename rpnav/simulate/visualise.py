from math import log2

import astropy.units as u
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from astropy.coordinates import ICRS, AltAz, SkyCoord
from astropy.table import Table
from astropy.time import Time

from rpnav.observe.antenna import Antenna
from rpnav.observe.pulsar import Pulsar

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
            marker = {'size': err.apply(lambda x: 1 / x)},
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
    fig.update_layout(
        xaxis_title="Longitude (deg)",
        yaxis_title="Longitude (deg)",
    )
    return fig



def plot_err_surface(x: list[float], y: list[float], residuals_rmse: list[float]):
    fig = go.Figure(data=[go.Surface(
        x=x,
        y=y,
        z=residuals_rmse,
        colorscale=[[0, 'rgb(0,100,80)'], 
              [1, 'rgba(68, 68, 68, 0.3)']]    
    )])
    fig.update_traces(
        contours_z=dict(show=True, usecolormap=True, highlightcolor="limegreen", project_z=True)
    )
    fig.update_layout(
        title="Timing Residuals RMSE Surface",
        # autosize=False,
        # scene_camera_eye=dict(x=1.87, y=0.88, z=-0.64),
        width=1800,
        height=1800,
        # margin=dict(l=65, r=50, b=65, t=90),
    )
    return fig