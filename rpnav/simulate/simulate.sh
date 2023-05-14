export TEMPO2=~/workspace/marsfield_sim/T2runtime/
ptaSimulate sim.input
tcsh msfd_sim/scripts/runScripts_master
tempo2 -output general2 -s "{sat} {pre} {err} rslt\n" -f msfd_sim/output/real_0/J0437-4715.par  msfd_sim/output/real_0/J0437-4715.tim | grep rslt | awk '{print $1,$2,$3}' > msfd_sim/output/real_0/results.txt
