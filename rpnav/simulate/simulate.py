import os
import subprocess
import io
from pathlib import Path
import shutil
import re

import numpy as np
import astropy.units as u

from rpnav.observe.antenna import Antenna
from rpnav.observe.observation import Observation
from rpnav.observe.pulsar import Pulsar
from rpnav.observe.observatories import PARKES
from rpnav.observe.observatories import WOODCHESTER
from astropy import constants as const

from rpnav.run import RunParams, runSimulate, run

FILE_DIR = Path(__file__).parent
INPUT_DIR = f"{FILE_DIR}/inputs/"
OUTPUT_DIR=f"{FILE_DIR}/outputs"

class SimParams(RunParams):
    def __init__(
        self, 
        name: str,
        antenna: Antenna,
        t: int,
        toa_err: dict[str: float],
        psr1: str,
        psr2: str,
        resultsDir: str,
        timFiles: dict[str, str] = {}, 
        parFiles: dict[str, str] = {},
        err: str = "chi",
        niter: int = 100,
        nsim: int = 100
    ):
        self.toa_err = toa_err
        super().__init__(
        name=name,
        antenna=antenna,
        t=t,
        psr1=psr1,
        psr2=psr2,
        resultsDir=resultsDir,
        timFiles=timFiles, 
        parFiles=parFiles,
        err=err,
        niter=niter,
        nsim=nsim
    )

def loadConfig(sectName: str, inputDict: dict[str, str], line: str) -> dict:
    """
    Load PTA Simulate Config
    """
    configName, configParams = line.split(":")
    if "=" in line:
        params = configParams.split(",")

        inputDict[sectName].append({configName: {}})

        for param in params:
            key, val = param.split("=")
            inputDict[sectName][-1][configName][key.strip()] = val.strip()
    else:
        inputDict[sectName].append({configName: configParams.strip()})
    return inputDict


def loadInput(inFile: str):
    """
    Load PTA Simulate file
    """
    inputDict = {}
    inSect = False
    sectName = ""
    with open(inFile, "r") as f:
        for line in f.readlines():
            line = line.strip()

            if len(line) == 0 or line.isspace() or line[0] == "#":
                continue
            elif line[0] == "<" and not inSect:
                sectName = line.replace('<', '').replace('>', '')
                inputDict[sectName] = []
                inSect = True
            elif line[0] == "<" and inSect:
                inSect = False
                sectName = ""
            elif inSect:
                inputDict = loadConfig(sectName, inputDict, line)
            else:
                exit("Error invalid input")
    return inputDict
            
def saveConfig(f: io.TextIOWrapper, config: dict[str, str]) -> dict:
    """
    Save PTA Simulate Config
    """
    for configName, config in config.items():
        if type(config) == dict:
            line = f"{configName}: "
            for k2, v2 in config.items():
                line += f"{k2}={v2},"
            # Remove trailing comma
            line = line[:-1]
            f.write(line)
        else:
            f.write(f"{configName}: {config}")

def saveInput(outputFile: str, inputDict):
    """
    Save PTA Simulate file
    """
    with open(outputFile, "w+") as f:
        for sectName, sectVal in inputDict.items():
            print(sectVal)
            f.write(f"<{sectName}>\n")
            if type(sectVal) == list:
                for config in sectVal:
                    saveConfig(f, config)
                    f.write(f"\n")
            f.write(f"</{sectName}>\n\n")
        
def formatInput(sim: SimParams, inputFile: str, outFile: str, simName: str):
    inputName = os.path.basename(inputFile)
    simInput = loadInput(inputFile)

    # Modify simulation name
    for i in range(len(simInput["define"])):
        defineDict = simInput["define"][i]
        if defineDict.get("name"):
            defineDict["name"]  = simName

    # Modify observed pulsars
    psrs = [sim.psr1, sim.psr2]
    for i in range(len(simInput["pulsars"])):
        defineDict = simInput["pulsars"][i]
        if defineDict.get("psr"):
            defineDict["psr"]["name"]  = psrs.pop()

    # Modify observation time
    startTime = None
    observer = ""
    for i in range(len(simInput["obsRun"])):
        defineDict = simInput["obsRun"][i]
        if defineDict.get("start"):
            startTime = float(defineDict["start"])
        elif defineDict.get("finish"):
            defineDict["finish"] = startTime + sim.t
        elif defineDict.get("name"):
            observer = defineDict.get("name")

    # Modify observation schedule
    psrs = [sim.psr1, sim.psr2]
    for i in range(len(simInput["schedule"])):
        defineDict = simInput["schedule"][i]
        if defineDict.get("observe"):
            psr = psrs.pop()
            defineDict["observe"]["psr"] = psr
            if sim.toa_err:
                defineDict["observe"]["toaerr"] = sim.toa_err[psr]

    print(simInput)
    saveInput(outFile, simInput)

def simulateTiming(sim: SimParams, outputFile: str, simName: str):
    runscript=f"{FILE_DIR}/{simName}/scripts/runScripts_master"

    ptaSimPath = os.environ.get("PTASIMULATE")
    tempo2Path = os.environ.get("TEMPO2")

    # Run PTA Simulate
    subprocess.run([
        f"{ptaSimPath} {outputFile}"
    ], shell=True)

    # Inject tempo2 path into run script
    subprocess.run([f"""echo 'alias tempo2 {tempo2Path}/bin/tempo2 
setenv TEMPO2 {tempo2Path}
' | \
cat - {runscript} >> temp \
&& mv temp {runscript}
"""
    ], shell=True)
    
    # Run simulate script
    subprocess.run([
        f"tcsh {runscript}"
    ], shell=True)
    shutil.move(f"{FILE_DIR}/{simName}", sim.resultsDir)


def run(sim: RunParams) -> dict:
    simName=f"{sim.name}_{sim.t}d_sim"
    simDir=f"{sim.resultsDir}/{simName}"
    simResDir=f"{simDir}/output/real_0"
    outdir=f"{sim.resultsDir}/{sim.psr1}/{sim.psr2}/{sim.t}h"
    os.makedirs(outdir, exist_ok=True)

    inputFile=f"{INPUT_DIR}/{sim.name}.input"
    outputFile = f"{outdir}/{simName}.input"

    # Clear old results
    print(f"{FILE_DIR}/{simName}")
    if os.path.exists(f"{FILE_DIR}/{simName}"):
        shutil.rmtree(f"{FILE_DIR}/{simName}")

    if os.path.exists(f"{sim.resultsDir}/{simName}"):
        shutil.rmtree(f"{sim.resultsDir}/{simName}")

    # Create input file
    formatInput(sim, inputFile, outputFile, simName)


    for si in range(1, sim.nsim):
        iterDir=f"{outdir}/S{si}"
        os.makedirs(iterDir, exist_ok=True)

        simulateTiming(sim, outputFile, simName)

        sim.parFiles[sim.psr1] = f"{simResDir}/{sim.psr1}.par"
        sim.timFiles[sim.psr1] = f"{simResDir}/{sim.psr1}.tim"
        sim.parFiles[sim.psr2] = f"{simResDir}/{sim.psr2}.par"
        sim.timFiles[sim.psr2] = f"{simResDir}/{sim.psr2}.tim"
        runSimulate(sim, iterDir) 
        shutil.rmtree(f"{sim.resultsDir}/{simName}")

# if __name__ == "__main__":
def test_main():
    observations = []

    # psrList=["J0835-4510", "J1017-7156", "J1024-0719", "J1600-3053", "J1732-5049", "J1909-3744", "J2129-5721", "J2241-5236", "J0613-0200", "J0711-6830", "J1022+1001", "J1045-4509", "J1125-6014", "J1446-4701", "J1545-4550", "J1603-7202", "J1643-1224", "J1713+0747", "J1730-2304", "J1744-1134", "J1824-2452A", "J1832-0836", "J1857+0943", "J1939+2134", "J2124-3358", "J2145-0750"]
    obsTime=[0.5, 1, 2, 3, 5, 8, 13, 21, 25, 28, 31]

    # wood_strong = Observation(
    #     observer=WOODCHESTER,
    #     pulsar=obs.pulsar,
    #     snr=20 * u.dimensionless_unscaled,
    #     integration_time=60 * 60 * u.s
    # )

    """
    Reference antenna defined in Feasibility Study for Pulsar Tim
    """
    SALA = Antenna(
        name="SALA",
        centre_frequency=1 * u.GHz,
        bandwidth=200 * u.MHz,
        effective_area=10 * u.m * u.m,
        location=PARKES.location,
        time=PARKES.time,
    )


    pulsars = Pulsar.load_catalogue()
    for p in pulsars:
        print(p.name)
        if p.name == "J1909-3744":
            observations.append(
                Observation(
                    PARKES,
                    p,
                    snr=146.09 * u.dimensionless_unscaled,
                    toa_err=(145 * u.ns).to(u.s),
                    integration_time= 1800 * u.s
                )
            )
        if p.name in ["J0835-4510", "B0833-45"]:
            observations.append(
                Observation(
                    PARKES,
                    p,
                    snr=146.09 * u.dimensionless_unscaled,
                    toa_err=(145 * u.ns).to(u.s),
                    integration_time= 1800 * u.s
                )
            )
        elif p.name == "B1937+21": 
            observations.append(
                Observation(
                    SALA,
                    p,
                    snr= (-55.6 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(3162 * u.m / const.c).to(u.s),
                    integration_time=(36.36 * u.min).to(u.s)
                )
            )
        elif p.name in ["B0736-40", "J0738-4042"]: 
            observations.append(
                Observation(
                    SALA,
                    p,
                    snr= (-50.2 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(9246467 * u.m / const.c).to(u.s),
                    integration_time=(3.07 * u.min).to(u.s)
                )
            )
        elif p.name == "B1451-68":
            observations.append(
                Observation(
                    SALA,
                    p,
                    snr= (-50.0 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(452529 * u.m / const.c).to(u.s),
                    integration_time=(2.35 * u.min).to(u.s)
                )
            )
        elif p.name == "B0950-08": 
            observations.append(
                Observation(
                    SALA,
                    p,
                    snr= (-49.8 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(316174 * u.m / const.c).to(u.s),
                    integration_time=(1.04 * u.min).to(u.s)
                )
            )
        elif p.name in ["B0329+54", "J0332+5434"]: 
            observations.append(
                Observation(
                    observer=SALA,
                    pulsar=p,
                    snr= (-45.2 * u.dB(u.dimensionless_unscaled)).to(u.dimensionless_unscaled),
                    toa_err=(290556 * u.m / const.c).to(u.s),
                    integration_time=(0.12 * u.min).to(u.s)
                )
            )

    psrNum = len(observations)

    observer: Antenna = WOODCHESTER
    print(observations)
    for t in obsTime:
        for pi in range(psrNum):
            psr1 = observations[pi].pulsar
            obs1 = Observation(
                observer=observer,
                pulsar=psr1,
                snr=5 * u.dimensionless_unscaled,
                integration_time=60 * 60 * u.s
            )

            for pj in range(psrNum):
                if pi >= pj:
                    continue    
                psr2 = observations[pj].pulsar

                obs2 = Observation(
                    observer=observer,
                    pulsar=psr2,
                    snr=5 * u.dimensionless_unscaled,
                    integration_time=60 * 60 * u.s
                )
                
                simName = f"{observer.name}_weak"
                RES_DIR=f"{OUTPUT_DIR}/{simName}"

                sim = SimParams(
                    name=simName,
                    antenna=observer,
                    psr1=psr1.name,
                    psr2=psr2.name,
                    resultsDir=RES_DIR,
                    toa_err={
                        psr1.name: float(obs1.toa_err(reference=obs1).value),
                        psr2.name: float(obs2.toa_err(reference=obs2).value)
                    },
                    t=t
                )
                run(sim)
                # run(sim)
                return

