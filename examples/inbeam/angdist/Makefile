all: s44_1329_histos.root   s44_1329_iso_histos.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina $*.mac | tee $*.log

%_histos.root : %.dat
	grutinizer gvalues.val -qH $*.dat -o /dev/null --hist-output $*_histos.root | tee $*_grut.log
	rm $*.dat

clean:
	rm *.log *.dat *.root
