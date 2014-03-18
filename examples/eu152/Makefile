all:  eu152_sim_histos.root eu152_gammas_sim_histos.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina $*.mac | tee $*.log

%_sim.root : %.dat
	SimCalculate -i $*.dat -o $*_sim.root -s eu152.set | tee $*_sim.log

#%_sim_histos.root : %_sim.root
#	Sim_histos -i $*_sim.root -o $*_sim_histos.root | tee $*_sim_histos.log

eu152_sim_histos.root : sm152_sim.root gd152_sim.root
	Sim_histos -i sm152_sim.root gd152_sim.root -o eu152_sim_histos.root | tee eu152_sim_histos.log

eu152_gammas_sim_histos.root : eu152_gammas_sim.root
	Sim_histos -i eu152_gammas_sim.root -o eu152_gammas_sim_histos.root | tee eu152_gammas_sim_histos.log

clean:
	rm *.log *.dat *.root