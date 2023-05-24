#Pulsar Timing Array (PTA) Simulation
TODO outline simulation paramaters 
## Usage
To set up Tempo2 local environment, 
```
export TEMPO2=/path/to/T2runtime
```

To set up the simulation environment, modify the `sim.input` file and run
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


