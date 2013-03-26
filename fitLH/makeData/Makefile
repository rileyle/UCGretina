all: example_ions.log \
     example_500.root \
     example_1000.root \
     example_1500.root \
     example_2000.root \
     example_spectrum.rad \
     example_spectrum.asc
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

example_spectrum.root : example_500.root   example_1000.root \
                        example_1500.root  example_2000.root
	root -l -b -q 'makeSpec.C'

example_spectrum.rad : example_spectrum.root
	root -l -b -q 'root2rad.C("example_spectrum.root", "Example Spectrum")'

example_spectrum.asc : example_spectrum.rad
	rad2asc2 example_spectrum.rad example_spectrum.asc

clean:
	rm *.log *.out *.events *.root *.spe *.wrl