# Pulsar Timing Array (PTA) Simulation
This package was designed to simulate timing observations received from radio pulsars. The

TODO outline simulation paramaters 
## Install
PTASimulate can be downloaded from [this repo](https://bitbucket.csiro.au/projects/PSRSOFT/repos/ptasimulate/). Once installed, modify the following parameter in make file to reference the ptasimulate location.
```
PREFIX := /path/to/ptasimulate
```
In the base directory run the following:
```
// build package
make

//install package 
make install
```

To make the executable available throughout the system, add the following line to `~/.bashrc`.
```
alias ptasimulate="/path/to/ptasimulate/ptaSimulate"
```

## Usage
To set up TEMPO2 local environment, 
```
export TEMPO2=/path/to/T2runtime
```

If the environment is set up according to the instructions in Dependencies.md, this should be `/usr/local/tempo2`.

## Simulate
To set up the simulation environment, configure the `sim.input` file to define the simulation scenario. A sample `.input` file is shown below
```
# Simulate Parkes telescope observing Vela 
<define>
name: parkes_sim            // Name of simulation
nproc: 1                    // Number of threads
nreal: 1
</define>

// Pulsar IDs 
<pulsars> 
psr: name=J0835-4510         
psr: name=J1022+1001        
psr: name=J1909-3744
</pulsars>

<obsRun>
name: PARKES                // Antenna name
tel: pks                    // Observatory ID defined in /simulate/T2runtime/observatory/
start: 60002.3              // Start JD time 
finish: 60003.3
sampling: cadence=0.1
sched: sched1
</obsRun>

<schedule>
name: sched1
observe: psr=J0835-4510,toaerr=1e-6,freq=1400,ha=8*(ran(linear)-0.5)
observe: psr=J0835-4510,toaerr=1e-6,freq=1400,ha=8*(ran(linear)-0.5)
observe: psr=J0835-4510,toaerr=1e-6,freq=1400,ha=8*(ran(linear)-0.5)
</schedule>
```


 and run the follow
```
ptaSimulate sim.input
```

```
# Run simulation
tcsh <script name>/scripts/runScripts_master

# output to results file
outpath=<script name>/output/real_<n>
tempo2 -output general2 -s "{sat} {pre} {err} rslt\n" -f ${outpath}/J0437-4715.par ${outpath}/J0437-4715.tim  | grep rslt | awk '{print $1,$2,$3}' > ${outpath}/results.txt

# convert from TCB to TDB
tempo2 -gr transform ${outpath}/J0835-4510.par ${outpath}/J0835-4510.tdb.par back
```


