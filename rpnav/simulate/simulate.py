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
from rpnav.observe.reference import OBSERVATIONS
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
            for key, param in config.items():
                line += f"{key}={param},"
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
            f.write(f"<{sectName}>\n")
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

    # Clear old simulation results
    if os.path.exists(f"{FILE_DIR}/{simName}"):
        shutil.rmtree(f"{FILE_DIR}/{simName}")

    if os.path.exists(f"{sim.resultsDir}/{simName}"):
        shutil.rmtree(f"{sim.resultsDir}/{simName}")

    # Create input file based on simulation parameters
    formatInput(sim, inputFile, outputFile, simName)

    for si in range(0, sim.nsim):
        iterDir=f"{outdir}/S{si}"
        os.makedirs(iterDir, exist_ok=True)

        simulateTiming(sim, outputFile, simName)

        sim.parFiles = {
            sim.psr1: f"{simResDir}/{sim.psr1}.par",
            sim.psr2: f"{simResDir}/{sim.psr2}.par",
        }
        sim.timFiles = {
            sim.psr1: f"{simResDir}/{sim.psr1}.tim",
            sim.psr2: f"{simResDir}/{sim.psr2}.tim",
        }
        runSimulate(sim, iterDir) 
        if os.path.exists(f"{sim.resultsDir}/{simName}"):
            shutil.rmtree(f"{sim.resultsDir}/{simName}")


if __name__ == "__main__":
    observer: Antenna = WOODCHESTER
    observations: list[Observation] = OBSERVATIONS
    obsTime=[1, 2, 3, 5, 8, 13, 21, 25, 28, 31]

    simName = f"{observer.name}_weak"
    RES_DIR=f"{OUTPUT_DIR}/{simName}"

    toa_errtol =  10 # 1e-3  
    psrs = list({obs.pulsar for obs in observations})
    psrNum = len(psrs)
    obsTime = [1]
    for t in obsTime:
        for pi, psr1 in enumerate(psrs):
            ref_obs1 = max([obs for obs in observations if obs.pulsar.name == psr1.name], key=lambda x: x.toa_err())

            obs1 = Observation(
                name=f"{psr1.name}_{simName}",
                observer=observer,
                pulsar=psr1,
                snr=5 * u.dimensionless_unscaled,
                integration_time=60 * 60 * u.s
            )
            obs1_toa_err = obs1.toa_err(reference=ref_obs1).to_value(u.s)

            if obs1_toa_err > toa_errtol:
                continue

            print(f"{ref_obs1.pulsar.name} & {ref_obs1.pulsar.jname} & {obs1_toa_err} \\\ ")
            continue
            for pj, psr2 in enumerate(psrs):
                if pi >= pj:
                    continue    

                ref_obs2 = max([obs for obs in observations if obs.pulsar.name == psr2.name], key=lambda x: x.toa_err())

                obs2 = Observation(
                    name=f"{psr2.name}_{simName}",
                    observer=observer,
                    pulsar=psr2,
                    snr=5 * u.dimensionless_unscaled,
                    integration_time=60 * 60 * u.s
                )
                obs2_toa_err = obs2.toa_err(reference=ref_obs2).to_value(u.s)

                if obs2_toa_err > toa_errtol:
                    continue

                sim = SimParams(
                    name=simName,
                    antenna=observer,
                    psr1=psr1.name,
                    psr2=psr2.name,
                    resultsDir=RES_DIR,
                    toa_err={
                        psr1.name: float(obs1.toa_err(reference=ref_obs1).to_value(u.s)),
                        psr2.name: float(obs2.toa_err(reference=ref_obs2).to_value(u.s))
                    },
                    t=t
                )
                run(sim)

