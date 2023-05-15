export TEMPO2=~/workspace/rpnav/rpnav/simulate/T2runtime
ptaSimulate sim.input
tcsh msfd_sim/scripts/runScripts_master
outpath=msfd_sim/output/real_0
tempo2 -output general2 -s "{sat} {pre} {err} rslt\n" -f ${outpath}/J0437-4715.par  ${outpath}/J0437-4715.tim | grep rslt | awk '{print $1,$2,$3}' > ${outpath}/results.txt
