all: s44_1329_histos.root \
     s44_2150_histos.root \
     s44_2457_histos.root
	$(MAKE) -C makeData
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina_LH $*.mac 2>&1 | tee $*.log
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
	$(MAKE) -C makeData clean

