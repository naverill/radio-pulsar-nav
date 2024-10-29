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
from rpnav.observe.observatories import PARKES, MSFD, WOODCHESTER_STRONG, WOODCHESTER_WEAK, NAVIGATE
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
        antenna = WOODCHESTER_WEAK
    )

# wood_strong_sim = "woodchester_weak"
wood_strong_sim = "woodchester_strong"
wood_strong = SimulationParams(
        name=wood_strong_sim,
        timfile=f"{FILE_DIR}/outputs/{wood_strong_sim}_sim/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/outputs/{wood_strong_sim}_sim/output/real_0/J0835-4510.tdb.par",
        antenna = WOODCHESTER_STRONG
    )

nav_sim = "navigate_small"
navigate = SimulationParams(
        name=nav_sim,
        timfile=f"{FILE_DIR}/outputs/{nav_sim}_sim/output/real_0/J0835-4510.tim",
        parfile=f"{FILE_DIR}/outputs/{nav_sim}_sim/output/real_0/J0835-4510.tdb.par",
        antenna = NAVIGATE
    )

@pytest.fixture(scope="module")
def sim() -> SimulationParams:
    return navigate


def fit_results(latlong, k: int = 2) -> KMeans:
    kmeans = KMeans(n_clusters=k)
    kmeans.fit(latlong)
    return kmeans


def test_rmse_surface(sim):
    fpath = f"{FILE_DIR}/outputs/parkes_1d_sim/output/"

    df1 = pd.read_csv(fpath + "ALPSMLC30_S033E148_DSM_GRID_CHI.csv", index_col=0)
    print(df1.shape)
    df1 = df1.iloc[:, [i for i in range(3600) if (i % 3 == 0)]]
    
    df2 = pd.read_csv(fpath + "ALPSMLC30_S034E148_DSM_GRID_CHI.csv", index_col=0)
    df2 = df2.iloc[:, [i for i in range(3600) if (i % 3 == 0)]]

    df = pd.concat([df1, df2])

    print(df.columns.values)
    df.columns.values.astype(float)
    long = df.columns.values.astype(float)
    lat = df.index.values
    err = df.values

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

def test_grde(sim):
    for i in range(100):
        simdir = f"S{i:02d}"
        fpath = f"{FILE_DIR}/outputs/{sim.name}/output/{simdir}/"

        # fname = "ALPSMLC30_N033E148_DSM_GRDE_CHI"

        fname = f"results_{simdir}"
        resname = f"{sim}_results_{simdir}"

        df = pd.read_csv(fpath + fname + ".csv")
        long = df["Longitude(deg)"].values
        lat = df["Latitude(deg)"].values

        fig= plot_scattermap(
                long, 
                lat, 
                df["Error"].values)
        fig.add_trace(go.Scatter(x=[148.26356], y=[-32.99841], marker_symbol='x', marker_size=40))

        fig.write_image(f"outputs/{resname}.png")

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
            print("ERR 1")
            continue

        for i in range(100):
            simdir = f"S{i:02d}"
            fpath = f"{simpath}/output/{simdir}/"
            if not os.path.exists(fpath):
                print("ERR 2")
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
    for t in [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]:
        print(t)
        sim_name = sim.name + f"_{t}d_sim"
    
        fig = go.Figure()
        fig.add_trace(
            go.Scattermap(
                lon=[float(sim.antenna.location.lon.to_value(u.deg))],
                lat=[float(sim.antenna.location.lat.to_value(u.deg))],
                marker={"color": "black", "size": 10, "symbol":"circle-x"},
                opacity=0.5,
            )
        )

        simpath = f"/media/shared/outputs/sim_data/{sim_name}"
        # fpath = f"{FILE_DIR}/outputs/{sim_name}/"
        if not os.path.exists(simpath):
            continue

        k = 2
        for i in range(100):
            simdir = f"S{i:02d}"
            fpath = f"{simpath}/output/{simdir}/"
            if not os.path.exists(fpath):
                continue

            fname = f"results_{simdir}"

            df = pd.read_csv(fpath + fname + ".csv")
            print(i)
            kmodel = fit_results(df[["Longitude(deg)", "Latitude(deg)"]], k=k)

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
                        opacity=0.2,
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

        # Create and add slider
        steps = []
        for i in range(len(fig.data)):
            step = dict(
                method="update",
                args=[{"visible": [False] * len(fig.data)},
                    {"title": "Slider switched to step: " + str(i)}],  # layout attribute
            )
            step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
            steps.append(step)


        # fig.write_image(f"outputs/{resname}_zoom.png")
        # fig.update_layout(map = {'zoom': 2.75, 'center': {'lon': 0, 'lat': 0}})
        # fig.write_image(f"outputs/{resname}.png")
        fig.show()
        # return



def test_init_scattermap(sim):
    sim_name = sim.name
    fig = go.Figure()
    fig.add_trace(
        go.Scattermap(
            lon=[float(sim.antenna.location.lon.to_value(u.deg))],
            lat=[float(sim.antenna.location.lat.to_value(u.deg))],
            marker_symbol='x', 
            marker={"color": "black", "size": 40}
        )
    )

    for i in range(100):
        simdir = f"S{i:02d}"
        lat = []
        long = []
        err = []
        for i in range(100):
            simfile = f"{simdir}I{i:02d}"
            fname = f"results_{simfile}"
            fpath = f"{FILE_DIR}/outputs/{sim_name}/output/{simdir}/{fname}.csv"
            if not os.path.exists(fpath):
                print("ERROR")
                continue

            df = pd.read_csv(fpath)
            long.append(df["Longitude(deg)"].values[0])
            lat.append(df["Latitude(deg)"].values[0])
            err.append(df["Error"][0])

        fig.add_trace(
            go.Scattermap(
                mode = "markers",
                marker = {'size': 10},
                lon = long,
                lat = lat,
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

def test_time_scattermap(sim):
    fig = go.Figure()
    fig.add_trace(
        go.Scattermap(
            lon=[float(sim.antenna.location.lon.to_value(u.deg))],
            lat=[float(sim.antenna.location.lat.to_value(u.deg))],
            marker_symbol='x', 
            marker={"color": "black", "size": 40}
        )
    )

    for t in [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]:
        lat = []
        long = []
        err = []
        sim_name = sim.name + f"_{t}d_sim"

        for i in range(100):
            simdir = f"S{i:02d}"
            fpath = f"{FILE_DIR}/outputs/{sim_name}_sim/output/{simdir}/"
            if not os.path.exists(fpath):
                continue

            fname = f"results_{simdir}"
            resname = f"{sim_name}_results_{simdir}"

            df = pd.read_csv(fpath + fname + ".csv")
            long.extend(list(df["Longitude(deg)"].values))
            lat.extend(list(df["Latitude(deg)"].values))
            err.extend(list(df["Error"]))

        fig.add_trace(
            go.Scattermap(
                mode = "markers",
                marker = {'size': [(1 / x) for x in err]},
                lon = long,
                lat = lat,
                line={'width':1},
                opacity=1,
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
    fig = make_subplots(rows=2, cols=2)
    for sim in [parkes, wood_weak, wood_strong]:
        XList = []
        YList = []
        ZList = []
        XStdList = []
        YStdList = []
        ZStdList = []
        TDList = []
        tlist = [0.168,0.33,0.5,1,2,3,5,8,13,21,25,28,31]
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

            k = 2
            for i in range(100):
                resid = 0
                simdir = f"S{i:02d}"
                fpath = f"{simpath}/output/{simdir}/"
                if not os.path.exists(fpath):
                    continue

                fname = f"results_{simdir}"

                df = pd.read_csv(fpath + fname + ".csv")

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
                    # X_stdev = cX.nanstd()
                    # Y_stdev = cY.nanstd()
                    # Z_stdev = cZ.nanstd()

                    if any(np.isnan([x, y, z])):
                        continue
                    if z > 0:
                        continue
                                
                    DX.append((x - true_x.to_value(u.m))**2) 
                    DY.append((y - true_y.to_value(u.m))**2)
                    DZ.append((z - true_z.to_value(u.m))**2)

            if len(DX) <= 0:
                continue
            
            EX = np.sqrt(sum(DX) / len(DX))
            EY = np.sqrt(sum(DY) / len(DY))
            EZ = np.sqrt(sum(DZ) / len(DZ))

            XList.append(EX)
            YList.append(EX)
            ZList.append(EX)
            XStdList.append(np.nanstd(DX)) 
            YStdList.append(np.nanstd(DY))
            ZStdList.append(np.nanstd(DZ))
            TD = (EX + EY + EX) / 3
            print(sim_name)
            print(f"X (RMSE): {EX * u.m}")
            print(f"Y (RMSE): {EY * u.m}")
            print(f"Z (RMSE): {EZ * u.m}")
            print(f"Total Distance: {TD * u.m}")
            TDList.append(TD)


        fig.add_trace(go.Scatter(
            x=tlist,
            y=TDList,
            name=sim.name,
            ),
            row=1, col=1
        )
        fig.add_trace(go.Scatter(
            x=tlist,
            y=XList,
            name=sim.name,
            error_y=dict(
                type='data', # value of error bar given in data coordinates
                array=XStdList,
                visible=True)
            ),
            row=1, col=2
        )
        fig.add_trace(go.Scatter(
            x=tlist,
            y=YList,name=sim.name,
            error_y=dict(
                type='data', # value of error bar given in data coordinates
                array=YStdList,
                visible=True)
            ),
            row=2, col=1
        )
        fig.add_trace(go.Scatter(
            x=tlist,
            y=YList,name=sim.name,
            error_y=dict(
                type='data', # value of error bar given in data coordinates
                array=ZStdList,
                visible=True)
            ),
            row=2, col=2
        )
    fig.update_yaxes(type="log")
    fig.show()