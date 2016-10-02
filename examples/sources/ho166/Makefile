all:  ho166_histos.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina $*.mac | tee $*.log

%_histos.root : %.dat
	grutinizer -qH $*.dat -o $*.root --hist-output $*_histos.root
	rm $*.dat
	rm $*.root

clean:
	rm *.log *.dat *.root
