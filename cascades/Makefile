all:    parallel.root superDef.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.log : %.mac
	time UCGretina $*.mac > $*.log
	rm $*.out

%.out : %.mac
	if test -e $*.out; \
	then echo "Removing existing file $*.out"; rm $*.out; \
	fi
	time UCGretina $*.mac > $*.log

%.root : %.out
	root -l -b -q 'UCGretinaSort.C("$*.out")'

clean:
	rm *.log *.out *.root *.spe *.wrl
