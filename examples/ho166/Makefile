all:  ho166_sim_histos.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina $*.mac | tee $*.log

%_sim.root : %.dat
	SimCalculate -i $*.dat -o $*_sim.root -s ho166.set | tee $*_sim.log


%_sim_histos.root : %_sim.root
	Sim_histos -i $*_sim.root -o $*_sim_histos.root | tee $*_sim_histos.log

clean:
	rm *.log *.dat *.root