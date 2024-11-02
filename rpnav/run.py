import os
import subprocess
import io
from pathlib import Path

import numpy as np
import pandas as pd

from rpnav.observe.antenna import Antenna
from rpnav.observe.observatories import PARKES
from sklearn.cluster import KMeans

FILE_DIR = Path(__file__).parent
INPUT_DIR = f"{FILE_DIR}/simulate/inputs/ppta_data/averaged"
OUTPUT_DIR=f"{FILE_DIR}/simulate/outputs"
RES_DIR=f"{OUTPUT_DIR}/parkes_real"
TRUNC_DIR = f"{INPUT_DIR}/trunc"
K = 2

class RunParams:
    def __init__(
        self, 
        name: str,
        antenna: Antenna,
        t: int,
        psr1: str,
        psr2: str,
        resultsDir: str,
        timFiles: dict[str, str] = {}, 
        parFiles: dict[str, str] = {},
        err: str = "chi",
        niter: int = 100,
        nsim: int = 100
    ):
        self.name :str = name
        self.psr1 : str = psr1
        self.psr2 : str = psr2
        self.timFiles : dict[str, str] = parFiles
        self.parFiles : dict[str, str] = timFiles
        self.resultsDir : str = resultsDir
        self.antenna :Antenna = antenna
        self.err :str = err
        self.t : int = t
        self.niter : int = niter
        self.nsim : int = nsim

def truncateFile(sim: RunParams, inFile: str, outFile: str, psr: str):
    f = open(inFile, 'r')
    psrLines = f.readlines()
    psrLen = len(psrLines)
    if sim.t > psrLen:
        exit(f"ERROR: Invalid observation time {t} > {psrLen}")
    
    # Skip first two lines
    psrStart=np.random.randint(2, psrLen) 
    psrEnd=psrStart + sim.t

    while psrEnd > psrLen:
        print(f"ERROR: {psr} Invalid end time")
        psrStart=np.random.randint(2, psrLen) 
        psrEnd=psrStart + sim.t       

    with open(outFile, "w+") as fTrunc:
        fTrunc.write(psrLines[0])
        fTrunc.write(psrLines[1])

        for i in range(psrStart, psrEnd):
            fTrunc.write(psrLines[i])

    f.close()

def runIteration(sim: RunParams, simfile: int, resfile: str):
    tempo2Path = os.environ.get("TEMPO2")
    subprocess.run([
        f"{tempo2Path}/bin/tempo2 -gr pulsar_positioning -o {sim.antenna.name} "
        f"-f {sim.parFiles[sim.psr1]} {sim.timFiles[sim.psr1]} "
        f"-f {sim.parFiles[sim.psr2]} {sim.timFiles[sim.psr2]} " 
        f"-a grde -e {str(sim.err)} -n -r -s {simfile}"
    ], shell=True)

    with open(simfile, "r+") as fsim:
        with open(resfile, "a+") as fres:
            psrLines = fsim.readlines()

            fres.write(psrLines[-1])

def runSimulate(sim: RunParams, outdir: str):
    resultsFile=f"{outdir}/results_{sim.psr1}_{sim.psr2}_{sim.err}.csv"

    with open(resultsFile, "w+") as f:
        f.write("Iteration,Longitude(deg),Latitude(deg),X(m),Y(m),Z(m),Error,Step Size")
    
    for i in range(sim.niter):
        simFile=f"{outdir}/results_{sim.psr1}_{sim.psr2}_I{i}_{sim.err}.csv"
        runIteration(sim, simFile, resultsFile)

    df = pd.read_csv(resultsFile)
    try:
        kmeans = KMeans(n_clusters=K)
        kmeans.fit(df[["Longitude(deg)", "Latitude(deg)"]])
    except ValueError:
        return 
    
    for i in range(K):
        cent = kmeans.cluster_centers_
        cLong = df["Longitude(deg)"][kmeans.labels_ == i]
        cLat = df["Latitude(deg)"][kmeans.labels_ == i]
        lon = cent[i][0]
        lat = cent[i][1]

def run(sim: RunParams, psr1Par: str, psr1Tim: str, psr2Par: str, psr2Tim: str):
    tim1Name = os.path.basename(psr1Tim)
    tim2Name = os.path.basename(psr2Tim)

    sim.parFiles[sim.psr1] = psr1Par
    sim.parFiles[sim.psr2] = psr2Par
    for si in range(2, sim.nsim):
        outdir=f"{sim.resultsDir}/{sim.psr1}/{sim.psr2}/{sim.t}h/S{si}"
        os.makedirs(outdir, exist_ok=True)

        psr1TruncFile=f"{outdir}/{tim1Name}"
        truncateFile(sim, psr1Tim, psr1TruncFile, sim.psr1)
        sim.timFiles[sim.psr1] = psr1TruncFile

        psr2TruncFile=f"{outdir}/{tim2Name}"
        truncateFile(sim, psr2Tim, psr2TruncFile, sim.psr2)
        sim.timFiles[sim.psr2] = psr2TruncFile

        os.makedirs(outdir, exist_ok=True)
        runSimulate(sim, outdir) 


if __name__ == "__main__":
    psrList=["J0613-0200", "J0711-6830", "J1024-0719",  "J1017-7156", "J1022+1001", "J1024-0719"]
    obsTime=[744]

    psrNum = len(psrList)

    observer: Antenna = PARKES
    for t in obsTime:
        for pi in range(psrNum):
            psr1 = psrList[pi]

            for pj in range(psrNum):
                if pi >= pj:
                    continue    

                psr2 = psrList[pj]
                sim = RunParams(
                    name=observer.name + "_real",
                    antenna=observer,
                    psr1=psr1,
                    psr2=psr2,
                    resultsDir=RES_DIR,
                    t=t
                )
                psr1Par = f"{INPUT_DIR}/{psr1}.par"
                psr1Tim = f"{INPUT_DIR}/{psr1}.avg.tim"
                psr2Par = f"{INPUT_DIR}/{psr2}.par"
                psr2Tim = f"{INPUT_DIR}/{psr2}.avg.tim"

                run(sim, psr1Par, psr1Tim, psr2Par, psr2Tim)

