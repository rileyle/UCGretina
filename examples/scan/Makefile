all:  scan.out scan.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.out : %.mac
	UCGretina_Scan $*.mac 2>&1 | tee $*.log

%.root : %.out
	root -l -b -q 'ScanSort.C("$*.out")' 2>&1 | tee $*_sort.log

clean:
	rm *.log *.out *.root
