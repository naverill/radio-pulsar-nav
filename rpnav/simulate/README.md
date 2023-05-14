# Pulsar Timing Array (PTA) Simulation

## Usage

```
ptaSimulate <.dat file>

source <script name>/scripts/runScripts_master
cd <script name>/output/real_<n>/
tempo2 -gr plk -f J0437-4715.par J0437-4715.tim

# output to results file
tempo2 -output general2 -s "{sat} {pre} {err} rslt\n" -f J0437-4715.par J0437-4715.tim  | grep rslt | awk '{print $1,$2,$3}' > results.txt
```


