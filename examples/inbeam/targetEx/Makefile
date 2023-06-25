all: au197_548_histos.root    simple_279_histos.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina $*.mac 2>&1 | tee $*.log

simple_279_histos.root: simple_279.dat
	grutinizer -qH simple_279.dat -o /dev/null \
	${GRUTSYS}/lib/libMakeHistosSourcesSim.so \
	--hist-output simple_279_histos.root 2>&1 | tee simple_279_grut.log
	rm simple_279.dat

%_histos.root : %.dat
	grutinizer gvalues.val -qH $*.dat -o /dev/null \
	${GRUTSYS}/lib/libMakeHistosInBeam.so \
	--hist-output $*_histos.root 2>&1 | tee $*_grut.log
	rm $*.dat

clean:
	rm *.log *.root

