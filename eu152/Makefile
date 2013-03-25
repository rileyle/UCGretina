all:  eu152.out                     eu152.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations.
%.out : %.mac
	if test -e $*.out; \
	then echo "Removing existing file $*.out"; rm $*.out; \
	fi
	time UCGretina $*.mac > $*.log

# Pattern rule for sorting the simulation output into Root histograms
# and Radware ascii spectra
%.root : %.out	
	time root -l -b -q 'UCGretinaSourceSort.C("$*.out")'

clean:
	rm *.log *.out *.root *.rad *.spe *.wrl