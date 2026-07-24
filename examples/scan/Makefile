all:  scan.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.out : %.mac
	UCGretina_Scan $*.mac 2>&1 | tee $*.log
ifeq ($(G4MULTITHREADED), 1)
	./concat_out.sh $* 2>&1   | tee -a $*.log
endif

%.root : %.out
	root -l -b -q 'ScanSort.C("$*.out")' 2>&1 | tee $*_sort.log

clean:
	rm -f *.log *.out *.root
