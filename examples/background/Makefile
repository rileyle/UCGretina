all:  background_histos.root \
      muons_histos.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina $*.mac 2>&1 | tee $*.log

%_histos.root : %.dat
	grutinizer -qH $*.dat -o /dev/null --hist-output $*_histos.root 2>&1 | tee $*_grut.log
	rm $*.dat

clean:
	rm *.log *.root
