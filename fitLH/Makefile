all: example_ions.log \
     example_500.out     example_500.root  \
     example_1000.out    example_1000.root  \
     example_1500.out    example_1500.root  \
     example_2000.out    example_2000.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.log : %.mac
	time UCGretina_LH $*.mac > $*.log
	rm $*.out

%.out : %.mac
	if test -e $*.out; \
	then echo "Removing existing file $*.out"; rm $*.out; \
	fi
	time UCGretina_LH $*.mac > $*.log

%.root : %.out
	root -l -b -q 'UCGretinaSort.C("$*.out")'

clean:
	rm *.log *.out *.root *.spe *.wrl
