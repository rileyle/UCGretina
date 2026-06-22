all: au197_548_histos.root    simple_279_histos.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina $*.mac 2>&1 | tee $*.log
ifeq ($(G4MULTITHREADED), 1)
	./concat.sh $* 2>&1   | tee -a $*.log
endif

%_histos.root : %.dat
	grutinizer gvalues.val -qH $*.dat -o /dev/null \
	${GRUTSYS}/lib/libMakeHistosInBeam.so \
	--hist-output $*_histos.root 2>&1 | tee $*_grut.log
	rm $*.dat

clean:
	rm -f *.log *.dat *.root

