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


parkes_sim = "parkes"
parkes = SimulationParams(
        name=parkes_sim,
        timfile=f"{FILE_DIR}/outputs/{parkes_sim}/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/outputs/{parkes_sim}/output/real_0/J0835-4510.tdb.par",
        antenna = PARKES
    )


msfd_sim = "msfd_est"
msfd =  SimulationParams(
        name=msfd_sim,
        timfile=f"{FILE_DIR}/outputs/{msfd_sim}/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/outputs/{msfd_sim}/output/real_0/J0835-4510.tdb.par",
        antenna = MSFD
    )

# nav_sim = "navigate_small_0.33d_sim"
wood_weak_sim = "woodchester_weak"
wood_weak = SimulationParams(
        name=wood_weak_sim,
        timfile=f"{FILE_DIR}/outputs/{wood_weak_sim}/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/outputs/{wood_weak_sim}/output/real_0/J0835-4510.tdb.par",
        antenna = WOODCHESTER
    )

# wood_strong_sim = "woodchester_weak"
wood_strong_sim = "woodchester_strong"
wood_strong = SimulationParams(
        name=wood_strong_sim,
        timfile=f"{FILE_DIR}/outputs/{wood_strong_sim}_sim/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/outputs/{wood_strong_sim}_sim/output/real_0/J0835-4510.tdb.par",
        antenna = WOODCHESTER
    )

nav_sim = "navigate_small"
navigate = SimulationParams(
        name=nav_sim,
        timfile=f"{FILE_DIR}/outputs/{nav_sim}_sim/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/outputs/{nav_sim}_sim/output/real_0/J0835-4510.tdb.par",
        antenna = NAVIGATE
    )

parkes_real_sim = "parkes_real"
parkes_real = SimulationParams(
        name=parkes_real_sim,
        timfile=f"{FILE_DIR}/outputs/{parkes_real_sim}_sim/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/outputs/{parkes_real_sim}_sim/output/real_0/J0835-4510.tdb.par",
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


def test_err_heatmap(sim):
    fpath = f"{FILE_DIR}/../observe/outputs/"
    df1 = pd.read_csv(fpath + "ALPSMLC30_S033E148_DSM_GRID_CHI.csv", index_col=0)
    df1 = df1.iloc[:, [i for i in range(3600) if (i % 3 == 0)]]
    
    df2 = pd.read_csv(fpath + "ALPSMLC30_S034E148_DSM_GRID_CHI.csv", index_col=0)
    df2 = df2.iloc[:, [i for i in range(3600) if (i % 3 == 0)]]

    df = pd.concat([df1, df2])

    # df.columns.values = df.columns.values.to(float)
    # lat = df.index.values
    # err = df.values
    # print(lat)
    fig = plot_heatmap(df)
    fig.update_layout(
        title=dict(
            text=f"Reduced χ² Error",
            xanchor="center",
            yanchor="top",
            x=0.5,
        ),
    )
    
    fig.update_scenes(xaxis_title_text="Longitude (deg)",  
                      yaxis_title_text="Latitude (deg)",  
                      )
    # fig.write_image("outputs/residErrSurface.png")
    fig.show()

def test_err_surface(sim):
    fpath = f"{FILE_DIR}/outputs/parkes_1d_sim/output/"
    fpath = f"{FILE_DIR}/inputs/ppta_data/averaged/"
    df1 = pd.read_csv(fpath + "ALPSMLC30_S033E148_DSM_GRID_CHI.csv", index_col=0)
    df1 = df1.iloc[:, [i for i in range(3600) if (i % 3 == 0)]]
    
    df2 = pd.read_csv(fpath + "ALPSMLC30_S034E148_DSM_GRID_CHI.csv", index_col=0)
    df2 = df2.iloc[:, [i for i in range(3600) if (i % 3 == 0)]]

    df = pd.concat([df1, df2])

    long = df.columns.values.astype(float)
    lat = df.index.values
    err = df.values
    print(lat)
    fig = plot_err_surface(long, lat, err)
    fig.update_layout(
        title=dict(
            text=f"Reduced χ² Error Surface",
            xanchor="center",
            yanchor="top",
            x=0.5,
        ),
    )
    
    fig.update_scenes(xaxis_title_text="Longitude (deg)",  
                      yaxis_title_text="Latitude (deg)",  
                      zaxis_title_text="χ²")
    # fig.write_image("outputs/residErrSurface.png")
    fig.show()

def test_grde(parkes_data):
    psr_name = "J0613-0200_J1022+1001"
    t = "744h"
    simdir = f"S0"
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

        fpath = f"{FILE_DIR}/outputs/{parkes_data.name}/{psr_name}/{t}/{simdir}/results_{psr_name}_{t}_I{i}d_rms.csv"
        if not os.path.exists(fpath):
            continue

        df = pd.read_csv(fpath)
        long = df["Longitude(deg)"].values
        lat = df["Latitude(deg)"].values
        err = df["Error"]
        print(err)
        fig.add_trace(
            go.Scattermap(
                mode="lines",
                # marker = {'size': [5 * (1 / x) if x > 0.01 else 0.01 for x in err]},
                lon = long,
                lat = lat,
                # line={'width':1},
                # opacity=1,
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

def test_scattermap(sim):
    for t in [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]:
        print(t)
        sim_name = sim.name + f"_{t}d_sim"

        fig = go.Figure()
        fig.add_trace(
            go.Scattermap(
                lon=[float(sim.antenna.location.lon.to_value(u.deg))],
                lat=[float(sim.antenna.location.lat.to_value(u.deg))],
                marker_symbol='x', 
                marker={"color": "black", "size": 40}
            )
        )

        simpath = f"{FILE_DIR}/outputs/{sim_name}/"
        simpath = f"/media/shared/outputs/{sim_name}"
        if not os.path.exists(simpath):
            continue

        for i in range(100):
            simdir = f"S{i:02d}"
            fpath = f"{simpath}/output/{simdir}/"
            if not os.path.exists(fpath):
                continue

            fname = f"results_{simdir}"

            df = pd.read_csv(fpath + fname + ".csv")
            long = df["Longitude(deg)"].values
            lat = df["Latitude(deg)"].values
            err = df["Error"]

            fig.add_trace(
                go.Scattermap(
                    mode = "markers",
                    marker = {'size': [5 * (1 / x) if x > 0.01 else 0.01 for x in err]},
                    lon = long,
                    lat = lat,
                    line={'width':1},
                    opacity=1,
                )
            )

        fig.update_layout(
            margin ={'l':0,'t':0,'b':0,'r':0},
            width=3600,
            height=1800,
            map = {
                # 'style': "white-bg",
                'style': "open-street-map",
                'center': {'lon': 148.26356, 'lat': -32.99841},
                'zoom': 10})

        # fig.write_image(f"outputs/{resname}_zoom.png")
        # fig.update_layout(map = {'zoom': 2.75, 'center': {'lon': 0, 'lat': 0}})
        # fig.write_image(f"outputs/{resname}.png")
        fig.show()

def test_centroid_scattermap(sim):
    # for t in [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]:
    test = [
        ("J0613-0200", "J0711-6830"), 
        ("J0437-4715", "J0613-0200"), 
        ("J0437-4715", "J0711-6830"), 
        ("J1909-3744", "B0329+54"), 
        ("J1909-3744", "B1937+21"), 
        ("J1909-3744", "B1911+1347"), 
        ("B1911+1347", "B0329+54"), 
        ("B1911+1347", "B1937+21"), 
        ("B1911+1347", "B1738+0333"), 
        # ("J0613-0200", "J1022+1001"), 
        # ("J0613-0200", "J1024-0719"), 
        # ("J0711-6830", "J1024-0719"), 
    ]
    for t in [0.168, 1]:
        for pi, pj in test:
            print(t)

            fig = go.Figure()
            fig.add_trace(
                go.Scattermap(
                    lon=[float(sim.antenna.location.lon.to_value(u.deg))],
                    lat=[float(sim.antenna.location.lat.to_value(u.deg))],
                    marker={"color": "black", "size": 10, "symbol":"circle-x"},
                    opacity=0.5,
                )
            )

            # simpath = f"/media/shared/outputs/sim_data/{sim_name}"
            simpath = f"{FILE_DIR}/outputs/{sim.name}/{pi}/{pj}/{t}h"
            if not os.path.exists(simpath):
                continue

            k = 2
            for i in range(100):
                simdir = f"S{i:02d}"
                fpath = f"{simpath}/{simdir}/"
                if not os.path.exists(fpath):
                    print("ERROR", fpath)
                    continue

                fname = f"results_{pi}_{pj}_chi.csv"

                df = pd.read_csv(fpath + fname)
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
                                # 'color' : errScaled,
                                # 'colorbar_title':"Reduced Chi-squared error",
                                # 'colorscale': 'Blues',
                            },
                            lon = cLong,
                            lat = cLat,
                            opacity=0.5,
                        )
                    )


            fig.update_layout(
                # margin ={'l':0,'t':0,'b':0,'r':0},
                width=3600,
                height=1800,
                map = {
                    # 'style': "white-bg",
                    'style': "open-street-map",
                    'center': {'lon': 148.26356, 'lat': -32.99841},
                    'zoom': 10})

            # fig.write_image(f"outputs/{resname}_zoom.png")
            # fig.update_layout(map = {'zoom': 2.75, 'center': {'lon': 0, 'lat': 0}})
            # fig.write_image(f"outputs/{resname}.png")
            fig.show()

def test_centroid_real_scattermap(parkes_data):
    test = [
        ("J0613-0200", "J0711-6830"), 
        # ("J0613-0200", "J1022+1001"), 
        # ("J0613-0200", "J1024-0719"), 
        # ("J0711-6830", "J1024-0719"), 
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
        # simpath = f"/media/shared/outputs/sim_data/{sim_name}"
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
                            # 'color' : errScaled,
                            # 'colorbar_title':"Reduced Chi-squared error",
                            # 'colorscale': 'Blues',
                        },
                        lon = cLong,
                        lat = cLat,
                        opacity=0.2,
                    )
                )
        print("Rate:", lenf / lenUnf)


        fig.update_layout(
            # margin ={'l':0,'t':0,'b':0,'r':0},
            width=3600,
            height=1800,
            map = {
                # 'style': "white-bg",
                'style': "open-street-map",
                'center': {'lon': 148.26356, 'lat': -32.99841},
                'zoom': 10})

        # fig.write_image(f"outputs/{resname}_zoom.png")
        # fig.update_layout(map = {'zoom': 2.75, 'center': {'lon': 0, 'lat': 0}})
        # fig.write_image(f"outputs/{resname}.png")
        fig.show()
        # return

def test_distance_err():
    fig = go.Figure()
    for sim in [parkes, wood_weak, wood_strong]:
        errList = []
        tlist = [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]
        for t in tlist:
            sim_name = sim.name + f"_{t}d_sim"
            lat, lon, _ = loc = sim.antenna.location.to_geodetic()
            lon = lon.to_value(u.deg)
            lat = lat.to_value(u.deg)

            true_x, true_y, true_z = sim.antenna.location.geocentric
            DX = []
            DY = []
            DZ = []
            for i in range(100):
                resid = 0
                simdir = f"S{i:02d}"
                # fpath = f"{FILE_DIR}/outputs/{sim_name}/output/{simdir}/"
                fpath = f"{FILE_DIR}/demo/{sim_name}/output/{simdir}/"
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
    test = ["J0613-0200_J0711-6830", "J0613-0200_J1022+1001", "J0613-0200_J1024-0719"]
    test = [
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
    for pi, pj in test:

        DX = []
        DY = []
        DZ = []
        print(t)

        # simpath = f"/media/shared/outputs/sim_data/{sim_name}"
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