# Excited states of 152Sm populated in 152Eu beta+/EC
SM = sm152_122_histos.root \
     sm152_366_histos.root \
     sm152_810_histos.root \
     sm152_1023_histos.root \
     sm152_1041_histos.root \
     sm152_1086_histos.root \
     sm152_1234_histos.root \
     sm152_1293_histos.root \
     sm152_1372_histos.root \
     sm152_1530_histos.root \
     sm152_1579_histos.root \
     sm152_1613_histos.root \
     sm152_1650_histos.root \
     sm152_1730_histos.root \
     sm152_1757_histos.root \
     sm152_1769_histos.root \
     sm152_1779_histos.root

# Excited states of 152Gd populated in 152Eu beta+/EC
GD = gd152_344_histos.root \
     gd152_755_histos.root \
     gd152_931_histos.root \
     gd152_1109_histos.root \
     gd152_1123_histos.root \
     gd152_1282_histos.root \
     gd152_1318_histos.root \
     gd152_1434_histos.root \
     gd152_1550_histos.root \
     gd152_1606_histos.root \
     gd152_1643_histos.root \
     gd152_1692_histos.root

all:  eu152_gammas_histos.root \
      $(SM) $(GD) \
      eu152_histos.root
	echo "Geant4 simulations done."

# Pattern rules for running Geant4 simulations and sorting the output
# into Root histograms
%.dat : %.mac
	UCGretina $*.mac | tee $*.log

%_histos.root : %.dat
	grutinizer -qH $*.dat -o $*.root --hist-output $*_histos.root
	rm $*.dat
	rm $*.root

eu152_histos.root : $(SM) $(GD)
	gadd eu152_histos.root $(SM) $(GD)

clean:
	rm *.log *.dat *.root
