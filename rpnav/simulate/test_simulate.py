from pathlib import Path
import csv
import os

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

import pandas as pd
import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import EarthLocation, Longitude, Latitude, Angle
from astropy.time import Time
from sklearn.cluster import KMeans

from rpnav.simulate.visualise import plot_heatmap, plot_scattermap, plot_err_surface
from rpnav.observe.observer import Observer
from rpnav.observe.observatories import PARKES, MSFD, WOODCHESTER, NAVIGATE
from rpnav.observe.antenna import Antenna

FILE_DIR = Path(__file__).parent


class SimulationParams:
    def __init__(
        self, antenna: Antenna, timfile: str, parfile: str, name: str, 
    ):
        self.name :str = name
        self.timfile :str = timfile
        self.parfile :str = parfile
        self.antenna :Antenna = antenna

parkes = SimulationParams(
    name= "parkes",
    antenna = PARKES
)

wood_weak = SimulationParams(
    name="woodchester_weak",
    antenna = WOODCHESTER
)

wood_strong = SimulationParams(
    name= "woodchester_strong",
    antenna = WOODCHESTER
)

navigate = SimulationParams(
    name="navigate_small",
    antenna = NAVIGATE
)

parkes_real = SimulationParams(
    name="parkes_real",
    antenna = PARKES
)

@pytest.fixture(scope="module")
def sim() -> SimulationParams:
    return wood_weak

@pytest.fixture(scope="module")
def parkes_data() -> SimulationParams:
    return parkes_real


def fit_results(latlong, k: int = 2) -> KMeans:
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(latlong)
    return kmeans


def test_grde(parkes_data):
    psrs = [
        ("J0613-0200", "J1022+1001"), 
    ]
    t = "744h"
    simdir = f"S0"
    err = "chi"
    for pi, pj in psrs:
        for i in range(100):
            fig = go.Figure()
            fig.add_trace(
                go.Scattermap(
                    lon=[float(parkes_data.antenna.location.lon.to_value(u.deg))],
                    lat=[float(parkes_data.antenna.location.lat.to_value(u.deg))],
                    marker={"color": "black", "size": 10, "symbol":"circle-x"},
                    opacity=0.5,
                )
            )

            fpath = f"{FILE_DIR}/outputs/{parkes_data.name}/{pi}/{pj}/{t}/{simdir}/results_{pi}_{pj}_{t}_I{i}d_{err}.csv"
            if not os.path.exists(fpath):
                continue

            df = pd.read_csv(fpath)
            long = df["Longitude(deg)"].values
            lat = df["Latitude(deg)"].values
            err = df["Error"]

            fig.add_trace(
                go.Scattermap(
                    mode="lines",
                    lon = long,
                    lat = lat,
                    text=[str(e) for e in err] 
                )
            )

            fig.update_layout(
                margin ={'l':0,'t':0,'b':0,'r':0},
                width=3600,
                height=1800,
                map = {
                    'style': "open-street-map",
                })
            fig.show()

def test_centroid_scattermap(sim):
    psrs = [
        ("J0613-0200", "J0711-6830"), 
    ]
    k = 2
    for t in [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]:
        for pi, pj in psrs:
            fig = go.Figure()
            fig.add_trace(
                go.Scattermap(
                    lon=[float(sim.antenna.location.lon.to_value(u.deg))],
                    lat=[float(sim.antenna.location.lat.to_value(u.deg))],
                    marker={"color": "black", "size": 10, "symbol":"circle-x"},
                    opacity=0.5,
                )
            )

            simpath = f"{FILE_DIR}/outputs/{sim.name}/{pi}/{pj}/{t}h"
            if not os.path.exists(simpath):
                continue

            for i in range(100):
                simdir = f"S{i:02d}"
                fpath = f"{simpath}/{simdir}"
                if not os.path.exists(fpath):
                    continue

                fname = f"results_{pi}_{pj}_chi.csv"

                df = pd.read_csv(f"{fpath}/{fname}")
                try:
                    kmodel = fit_results(df[["Longitude(deg)", "Latitude(deg)"]], k=k)
                except ValueError:
                    continue

                for i in range(k):
                    cent = kmodel.cluster_centers_
                    cLong = df["Longitude(deg)"][kmodel.labels_ == i]
                    cLat = df["Latitude(deg)"][kmodel.labels_ == i]
                    lon = cent[i][0]
                    lat = cent[i][1]

                    fig.add_trace(
                        go.Scattermap(
                            lon=[lon],
                            lat=[lat],
                            marker={"color": "black", "size":5},
                            opacity=0.5,
                        )
                    )

                    err = df["Error"][kmodel.labels_ == i]

                    errScaled = [5 * (1 / x) if x > 0.01 else 0.01 for x in err]
                    fig.add_trace(
                        go.Scattermap(
                            mode = "markers",
                            marker = {
                                'size': errScaled,
                            },
                            lon = cLong,
                            lat = cLat,
                            opacity=0.5,
                        )
                    )

            fig.update_layout(
                margin ={'l':0,'t':0,'b':0,'r':0},
                width=3600,
                height=1800,
                map = {
                    'style': "open-street-map",
                    'center': {'lon': 148.26356, 'lat': -32.99841},
                    'zoom': 10})

            fig.show()

def test_centroid_real_scattermap(parkes_data):
    test = [
        ("J0613-0200", "J0711-6830"), 
    ]
    t = "744h"

    for pi, pj in test:
        print(t)
    
        fig = go.Figure()
        fig.add_trace(
            go.Scattermap(
                lon=[float(parkes_data.antenna.location.lon.to_value(u.deg))],
                lat=[float(parkes_data.antenna.location.lat.to_value(u.deg))],
                marker={"color": "black", "size": 10, "symbol":"circle-x"},
                opacity=0.5,
            )
        )

        lenUnf = 0
        lenf = 0
        simpath = f"{FILE_DIR}/outputs/parkes_real/{pi}/{pj}/{t}"
        if not os.path.exists(simpath):
            continue

        k = 2
        for i in range(100):
            simdir = f"S{i}"
            fpath = f"{simpath}/{simdir}/results_{pi}_{pj}_chi.csv"
            if not os.path.exists(fpath):
                print(fpath)
                print("ERROR")
                continue

            df = pd.read_csv(fpath)
            lenUnf += len(df["Error"])
            df = df[df["Error"] < 8000]
            lenf += len(df["Error"])
            try:
                kmodel = fit_results(df[["Longitude(deg)", "Latitude(deg)"]], k=k)
            except ValueError:
                continue

            for i in range(k):
                cent = kmodel.cluster_centers_
                cLong = df["Longitude(deg)"][kmodel.labels_ == i]
                cLat = df["Latitude(deg)"][kmodel.labels_ == i]
                lon = cent[i][0]
                lat = cent[i][1]

                fig.add_trace(
                    go.Scattermap(
                        lon=[lon],
                        lat=[lat],
                        marker={"color": 'red', "size":10},
                        opacity=0.5,
                    )
                )

                err = df["Error"][kmodel.labels_ == i]

                errScaled = [5 * (1 / x) if x > 0.01 else 0.01 for x in err]
                fig.add_trace(
                    go.Scattermap(
                        mode = "markers",
                        marker = {
                            'size': errScaled,
                        },
                        lon = cLong,
                        lat = cLat,
                        opacity=0.2,
                    )
                )
        print("Convergence Rate:", lenf / lenUnf)


        fig.update_layout(
            # margin ={'l':0,'t':0,'b':0,'r':0},
            width=3600,
            height=1800,
            map = {
                'style': "open-street-map",
                'center': {'lon': 148.26356, 'lat': -32.99841},
                'zoom': 10})

        fig.show()

def test_distance_err():
    fig = go.Figure()
    for sim in [parkes, wood_weak, wood_strong]:
        errList = []
        tlist = [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]
        for t in tlist:
            lat, lon, _ = loc = sim.antenna.location.to_geodetic()
            lon = lon.to_value(u.deg)
            lat = lat.to_value(u.deg)

            fpath = f"{FILE_DIR}/output/{sim.name}/output/{simdir}/"
            true_x, true_y, true_z = sim.antenna.location.geocentric
            DX = []
            DY = []
            DZ = []
            for i in range(100):
                resid = 0
                simdir = f"S{i:02d}"
                if not os.path.exists(fpath):
                    continue

                fname = f"results_{simdir}"
                df = pd.read_csv(fpath + fname + ".csv")

                if "X(m)" not in df.columns:
                    continue

                X = df["X(m)"].values
                Y = df["Y(m)"].values
                Z = df["Z(m)"].values
                lat = df["Latitude(deg)"].values
                err = df["Error"].values

                for i in range(len(err)):
                    x = X[i] * u.m
                    y = Y[i] * u.m
                    z = Z[i] * u.m
                    l = lat[i]
                    e = err[i]

                    if any(np.isnan([x.value, y.value, z.value])):
                        continue
                    if l > 0:
                        continue
                    if e > 2:
                        continue
                                
                    DX.append((x - true_x)**2) 
                    DY.append((y - true_y)**2)
                    DZ.append((z - true_z)**2) 

            if len(DX) <= 0:
                continue
            
            EX = np.sqrt(sum(DX) / len(DX))
            EY = np.sqrt(sum(DY) / len(DY))
            EZ = np.sqrt(sum(DZ) / len(DZ))
            TD = (EX + EY + EX) / 3
            print(sim_name)
            print(f"X (RMSE): {EX}")
            print(f"Y (RMSE): {EY}")
            print(f"Z (RMSE): {EZ}")
            print(f"Total Distance: {TD}")
            errList.append(TD.to_value(u.m))
        fig.add_trace(go.Scatter(x=tlist, y=errList, name=sim.name))
    fig.update_yaxes(type="log")
    fig.show()

def test_centroid_distance_err():
    fig = make_subplots(rows=3, cols=1)
    # for sim, lcol, fcol in [(wood_strong, 'rgba(0,0,0,0.5)','rgba(68, 68, 68, 0.3)'), (wood_weak, 'rgb(0,100,80)', 'rgba(0,100,80,0.2)')]:
    for sim, lcol, fcol in [(parkes, 'rgb(0,100,80)', 'rgba(0,100,80,0.2)')]:
        PXList = []
        PYList = []
        PZList = []
        NXList = []
        NYList = []
        NZList = []
        PXStdList = []
        PYStdList = []
        PZStdList = []
        NXStdList = []
        NYStdList = []
        NZStdList = []
        TDList = []
        tlist = [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]
        # tlist = [0.168,0.33,0.5,1,2,3,8,]
        for t in tlist:
            sim_name = sim.name + f"_{t}d_sim"

            # fpath = f"{FILE_DIR}/outputs/{sim_name}/"
            simpath = f"{FILE_DIR}/outputs/{sim_name}/"
            if not os.path.exists(simpath):
                continue

            true_x, true_y, true_z = sim.antenna.location.geocentric

            DX = []
            DY = []
            DZ = []
            TDX = []
            TDY = []
            TDZ = []
            k = 2
            for i in range(100):
                resid = 0
                simdir = f"S{i:02d}"
                fname = f"results_{simdir}"
                fpath = f"{simpath}/output/{simdir}/{fname}.csv"
                if not os.path.exists(fpath):
                    continue

                df = pd.read_csv(fpath)

                df = df[df["Error"] < 5] 

                if "X(m)" not in df.columns:
                    continue

                try:
                    kmodel = fit_results(df[["X(m)", "Y(m)", "Z(m)"]], k=k)
                except ValueError:
                    continue

                
                for c in range(k):
                    cent = kmodel.cluster_centers_
                    cX = df["X(m)"][kmodel.labels_ == c]
                    cY = df["Y(m)"][kmodel.labels_ == c]
                    cZ = df["Z(m)"][kmodel.labels_ == c]
                    x = cent[c][0]
                    y = cent[c][1]
                    z = cent[c][2]

                    if any(np.isnan([x, y, z])):
                        continue
                    if z > 0:
                        continue
                                
                    DX.append((x - true_x.to_value(u.m))) 
                    TDX.append(abs(x - true_x.to_value(u.m))) 
                    DY.append((y - true_y.to_value(u.m)))
                    TDY.append(abs(y - true_y.to_value(u.m)))
                    DZ.append((z - true_z.to_value(u.m)))
                    TDZ.append(abs(z - true_z.to_value(u.m)))

            if len(DX) <= 0:
                continue
            
            errX = np.nanmean(TDX)
            errY = np.nanmean(TDY)
            errZ = np.nanmean(TDZ)
            print()
            print(sim.name)
            print(t)
            print(f"Total X Error: {errX}")
            print(f"Total Y Error: {errY}")
            print(f"Total Z Error: {errZ}")

            PXList.append(np.nanmean([x for x in DX if x >= 0]))
            NXList.append(np.nanmean([x for x in DX if x < 0]))
            PYList.append(np.nanmean([x for x in DY if x >= 0]))
            NYList.append(np.nanmean([x for x in DY if x < 0]))
            PZList.append(np.nanmean([x for x in DZ if x >= 0]))
            NZList.append(np.nanmean([x for x in DZ if x < 0]))

            PXStdList.append(np.nanstd([x for x in DX if x >= 0])) 
            NXStdList.append(np.nanstd([x for x in DX if x < 0]))
            PYStdList.append(np.nanstd([x for x in DY if x >= 0]))
            NYStdList.append(np.nanstd([x for x in DY if x < 0])) 
            PZStdList.append(np.nanstd([x for x in DZ if x >= 0]))
            NZStdList.append(np.nanstd([x for x in DZ if x > 0]))

        print(PXStdList)
        print(NXStdList)
        print(PYStdList)
        print(NYStdList)
        print(PZStdList)
        print(NZStdList)
        ############### X
        fig.add_trace(go.Scatter(
                x=tlist,
                y=PXList,
                name=f"{sim.name.replace('_', ' ')} +X Error",
                mode='lines',
                marker=dict(color=lcol),

            ),
            row=1, col=1
        )
        fig.add_trace(go.Scatter(
                x=tlist,
                y=NXList,
                mode='lines',
                marker=dict(color=lcol),
                name=f"{sim.name.replace('_', ' ')} -X Error"
            ),
            row=1, col=1
        )
        fig.add_trace(go.Scatter(
                x=tlist,
                y=[m + std for m, std in zip(PXList, PXStdList)],
                name=f"{sim.name.replace('_', ' ')} Upper Bound",
                mode='lines',
                marker=dict(color=lcol),
                line=dict(width=0),
                showlegend=False
            ),
            row=1, col=1
        )

        fig.add_trace(go.Scatter(
                name='Lower Bound',
                x=tlist,
                y=[m - std for m, std in zip(NXList, NXStdList)],
                marker=dict(color=lcol),
                line=dict(width=0),
                mode='lines',
                fillcolor=fcol,
                fill='tonexty',
                showlegend=False
            ),
            row=1, col=1
        )
        ############### Y
        fig.add_trace(go.Scatter(
                x=tlist,
                y=PYList,
                mode='lines',
                marker=dict(color=lcol),
                name=f"{sim.name.replace('_', ' ')} +Y Error"
            ),
            row=2, col=1
        )
        fig.add_trace(go.Scatter(
                x=tlist,
                y=NYList,
                mode='lines',
                marker=dict(color=lcol),
                name=f"{sim.name.replace('_', ' ')} -Y Error"
            ),
            row=2, col=1
        )
        fig.add_trace(go.Scatter(
                x=tlist,
                y=[m + std for m, std in zip(PYList, PYStdList)],
                name=f"{sim.name.replace('_', ' ')} Upper Bound",
                mode='lines',
                marker=dict(color=lcol),
                line=dict(width=0),
                showlegend=False
            ),
            row=2, col=1
        )
        fig.add_trace(go.Scatter(
                name='Lower Bound',
                x=tlist,
                y=[m - std for m,std in zip(NYList, NYStdList)],
                marker=dict(color=lcol),
                line=dict(width=0),
                mode='lines',
                fillcolor=fcol,
                fill='tonexty',
                showlegend=False
            ),
            row=2, col=1
        )
        ############### Z
        fig.add_trace(go.Scatter(
                x=tlist,
                y=PZList,
                mode='lines',
                marker=dict(color=lcol),
                name=f"{sim.name.replace('_', ' ')} +Z Error"
            ),
            row=3, col=1
        )
        fig.add_trace(go.Scatter(
                x=tlist,
                y=NZList,
                mode='lines',
                marker=dict(color=lcol),
                name=f"{sim.name.replace('_', ' ')} -Z Error"
            ),
            row=3, col=1
        )
        fig.add_trace(go.Scatter(
                x=tlist,
                y=[m + std for m, std in zip(PZList, PZStdList)],
                name=f"{sim.name} Upper Bound",
                mode='lines',
                marker=dict(color=lcol),
                line=dict(width=0),
                showlegend=False
            ),
            row=3, col=1
        )
        fig.add_trace(go.Scatter(
                name='Lower Bound',
                x=tlist,
                y=[m - std for m , std in zip(NZList, NZStdList)],
                marker=dict(color=lcol),
                line=dict(width=0),
                mode='lines',
                fillcolor=fcol,
                fill='tonexty',
                showlegend=False
            ),
            row=3, col=1
        )

        # Update xaxis properties
        fig.update_yaxes(title_text="X (m)", type="log", row=1, col=1)
        fig.update_yaxes(title_text="Y (m)", type="log", row=2, col=1)
        fig.update_yaxes(title_text="Z (m)", type="log", row=3, col=1)

        # Update xaxis properties
        fig.update_xaxes(title_text="Observation time (days)", row=1, col=1)
        fig.update_xaxes(title_text="Observation time (days)", row=2, col=1)
        fig.update_xaxes(title_text="Observation time (days)", row=3, col=1)


        fig.update_layout(
            title=dict(
                text=f"Observatory Estimation Error",
                xanchor="center",
                yanchor="top",
                x=0.5,
            )
        )
    fig.show()  


def test_centroid_real_distance_err(parkes_data):
    psrs = [
        ("J0613-0200", "J0711-6830"), 
        # ("J0613-0200", "J1022+1001"), 
        # ("J0613-0200", "J1024-0719"), 
        # ("J0711-6830", "J1024-0719"), 
    ]
    t = "744h"

    lat, lon, _ = loc = parkes_data.antenna.location.to_geodetic()
    lon = lon.to_value(u.deg)
    lat = lat.to_value(u.deg)

    true_x, true_y, true_z = parkes_data.antenna.location.geocentric
    for pi, pj in psrs:

        DX = []
        DY = []
        DZ = []
        print(t)

        simpath = f"{FILE_DIR}/outputs/parkes_real/{pi}/{pj}/{t}"
        if not os.path.exists(simpath):
            continue

        k = 2
        for i in range(100):
            simdir = f"S{i:02d}"
            fpath = f"{simpath}/{simdir}/results_{pi}_{pj}_chi.csv"
            if not os.path.exists(fpath):
                continue

            
            df = pd.read_csv(fpath)
            df = df[df["Error"] < 500]

            if "X(m)" not in df.columns:
                continue

            try:
                kmodel = fit_results(df[["X(m)", "Y(m)", "Z(m)"]], k=k)
            except ValueError:
                continue

                
            for i in range(k):
                cent = kmodel.cluster_centers_
                cX = df["X(m)"][kmodel.labels_ == i]
                cY = df["Y(m)"][kmodel.labels_ == i]
                cZ = df["Z(m)"][kmodel.labels_ == i]
                x = cent[i][0]
                y = cent[i][1]
                z = cent[i][2]

                if any(np.isnan([x, y, z])):
                    continue
                if z > 0:
                    continue
                DX.append(abs(x - true_x.to_value(u.m))) 
                DY.append(abs(y - true_y.to_value(u.m)))
                DZ.append(abs(z - true_z.to_value(u.m)))
        errX = np.nanmean(DX)
        errY = np.nanmean(DY)
        errZ = np.nanmean(DZ)
        print()
        print(parkes_data.name)
        print(f"Total X Error: {errX}")
        print(f"Total Y Error: {errY}")
        print(f"Total Z Error: {errZ}")