all: s44_1329_histos.root \
     s44_2150_histos.root \
     s44_2457_histos.root
	$(MAKE) -C makeData
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina_LH $*.mac | tee $*.log

%_histos.root : %.dat
	grutinizer gvalues.val -qH $*.dat -o /dev/null --hist-output $*_histos.root | tee $*_grut.log
	rm $*.dat

clean:
	$(MAKE) -C makeData clean
	rm *.log *.root

