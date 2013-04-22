all: example_ions.log \
     example_500_sim_histos.root \
     example_1000_sim_histos.root  \
     example_1500_sim_histos.root  \
     example_2000_sim_histos.root
	$(MAKE) -C makeData
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.log : %.mac
	UCGretina_LH $*.mac > $*.log
	rm $*.out

%.dat : %.mac
	UCGretina_LH $*.mac > $*.log

%_sim.root : %.dat
	SimCalculate -i $*.dat -o $*_sim.root -s example.set > $*_sim.log

%_sim_histos.root : %_sim.root
	Sim_histos -i $*_sim.root -o $*_sim_histos.root > $*_sim_histos.log

clean:
	rm *.log *.dat *.root
	$(MAKE) -C makeData clean
