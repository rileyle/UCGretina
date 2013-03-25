#include "Riostream.h"
#include "iostream"
#include "TRandom.h"
#include "TFile.h"

TH1F *sim;
TH1F *sim_addback;
TH1F *crystal[28];
TH1F *ring[8];

Float_t beta, zoffset;
Float_t sigmaPar0, sigmaPar1, sigmaPar2, sigmaPar3;
Float_t ecalPar0, ecalPar1;
Float_t threshPar0[28] = {0.0};
Float_t threshPar1[28] = {0.0};

// 152Eu and 56Co data
const Int_t nPeaks = 20;
const Int_t iEu152 = 0; // Array index where 152Eu data begins
const Int_t nEu152 = 9; // Number of 152Eu efficiencies
const Int_t iCo56 = 9;  // Array index where 56Co data begins
const Int_t nCo56 = 11; // Number of 56Co efficiencies
Int_t iFirst, iLast;

// IAEA data
Float_t Energy[nPeaks] = {121.782, 244.699, 344.281,  411.126,  443.965, 778.903, 964.055, 1112.087, 1408.022, 846.764, 1037.884, 1175.099, 1238.287, 1360.206, 1771.350, 2034.759, 2598.460, 3201.954, 3451.154, 3548.27};
Float_t Branch[nPeaks] = { 0.2837,  0.0753,  0.2657,  0.02238,  0.03125,  0.1297,  0.1463,   0.1354,  0.2085,  0.99933,   0.1413,  0.02239,   0.6607,  0.04256,   0.1549,  0.07771,   0.1696,   0.0313,   0.0093, 0.00178};

Float_t gammasPerDecay[3] = {1.537, 2.435, 10.};

Float_t Eff[nPeaks] = {0.0};
Float_t EffAB[nPeaks] = {0.0};

Float_t gammasEmitted = 1e7;

// Read beta, target offset, and resolution parameters
void loadGretina(TString fileName) {

  char line[1000];
  Float_t fBuffer[4];
  TString sBuffer[8];

  ifstream fp;
  fp.open(fileName);

  // Beta for Doppler reconstruction
  fp >> beta;
  fp.getline(line,1000);  // Advance to next line.

  // Target z offset for Doppler reconstruction
  fp >> zoffset;
  fp.getline(line,1000);  // Advance to next line.

  // Resolution parameters
  fp >> sigmaPar0 >> sigmaPar1 >> sigmaPar2 >> sigmaPar3;
  fp.getline(line,1000);  // Advance to next line.

  // Energy calibration parameters
  fp >> ecalPar0  >> ecalPar1;
  fp.getline(line,1000);  // Advance to next line.

  // Threshold parameters
  Int_t i = 0;
  while(fp >> threshPar0[i] >> threshPar1[i]){
    cout << "Crystal " << i << " threshold parameters: "
	 << threshPar0[i] << ", " << threshPar1[i] << endl;
    fp.getline(line,1000);  // Advance to next line.
    i++;
  }
  if(i == 0)
    cout << "No thresholds set!" << endl;

  fp.close();

}

// Simulate intrinsic detector energy resolution and energy (de)calibration
Float_t measuredE(Float_t E){

  E = E*ecalPar1 + ecalPar0;

  Float_t sigma = sigmaPar0 
    + sigmaPar1*E
    + sigmaPar2*E*E
    + sigmaPar3*E*E*E;

  return E + gRandom->Gaus(0,sigma);
}

// Read simulation data and sort into a Doppler-corrected histogram
void loadSim(TString fileName, TString sourceType) {

  if(sourceType.Contains("eu152")){
    iFirst = iEu152;
    iLast  = iEu152 + nEu152;
  } else if(sourceType.Contains("co56")){
    iFirst = iCo56;
    iLast  = iCo56 + nCo56;
  } else if(sourceType.Contains("cs137")||sourceType.Contains("co60")){
    iFirst = 0;
    iLast  = 0;
  } else if(sourceType.Contains("photopeaks")){
    iFirst = 0;
    iLast  = nPeaks;
  } else {
    cout << "Error: Unknown source type " << sourceType << endl;
    return;
  }

  char line[1000];
  Float_t fBuffer[4];
  Int_t   iBuffer;
  TString sBuffer;

  Float_t Edep[1000] = {0.0};
  Int_t   detNum[200] = {0};

  ifstream fp;
  fp.open(fileName);

  Int_t nChannels     = 8192;
  Int_t keVperChannel = 1;    // Match binning of the spectrum to fit
  Int_t lo = 0;
  Int_t hi = lo+nChannels*keVperChannel;

  TString Name = fileName.Copy();
  Name.Replace(Name.Index(".out",4,0,0),4,"",0);
  sim = new TH1F(Name,"", nChannels, float(lo), float(hi));

  TString NameAB = fileName.Copy();
  NameAB.Replace(NameAB.Index(".out",4,0,0),4,"_addback",8);
  sim_addback = new TH1F(NameAB,"", nChannels, float(lo), float(hi));

  Int_t nCrystals = 24;
  // Crystal ID's in the Berkeley scheme: Hole = (int)ID/4, Crystal = ID%4
  //  Int_t crystalID[] = {24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,60,61,62,63,68,69,70,71};
  //  Int_t crystalID[] = {20,21,22,23, 24,25,26,27, 28,29,30,31, 32,33,34,35, 56,57,58,59, 64,65,66,67};

  // To match SpecTcl
  Int_t crystalID[] = {22,20,21,23, 25,27,26,24, 29,31,30,28, 33,35,34,32, 57,59,58,56, 65,67,66,64};
  Int_t crystalNum[100];

  for(Int_t i = 0; i < nCrystals; i++){
    // Load detector ID lookup table
    crystalNum[crystalID[i]] = i; 

    // Initialize crystal spectra
    TString crystalName = Name.Copy();
    TString crystalLabel;
    crystalLabel.Form("_%02d", i);
    crystalName += crystalLabel;
    crystal[i] = new TH1F(crystalName,"", nChannels, float(lo), float(hi));
    crystal[i]->Sumw2();
  }

  // Initialize ring spectra
  Int_t nRings = 8;
  Int_t ringAngle[] = {49, 55, 64, 67, 74, 87, 93, 105};
  for(Int_t i = 0; i < nRings; i++){
    TString ringName = Name.Copy();
    TString ringLabel;
    ringLabel.Form("_%02d", ringAngle[i]);
    ringName += ringLabel;
    ring[i] = new TH1F(ringName,"", nChannels, float(lo), float(hi));
  }

  Int_t nGamma = 0;
  Int_t nReaction = 0;
  Int_t nPhotopeak = 0;
  Int_t nGamma, nEvent;
  Int_t nHits;
  Float_t Egamma;

  while( fp >> sBuffer >> nHits >> Egamma >> nEvent){

    //    cout << sBuffer << "   " 
    //         << nHits << "   " 
    //         << Egamma << "   "
    //         << nEvent
    //         << endl;

    Float_t ElabAB = 0;
    Float_t x,   y, z;
    Float_t xs, ys, zs;
    Float_t EdepMax = 0;

    Float_t Edeposited;

    if(nHits > 200)
      cout << "\nError: nHits = " << nHits << endl;

    // Load the Edep array: energy deposited, indexed by detector number
    //    and detNum array: detector number, indexed by hit number
    for(Int_t i = 0; i < nHits; i++){

      fp >> iBuffer
	 >> fBuffer[0] >> fBuffer[1] >> fBuffer[2] >> fBuffer[3];

      //      cout << iBuffer 
      //      	   << "   " << fBuffer[0]
      //      	   << "   " << fBuffer[1]
      //      	   << "   " << fBuffer[2]
      //      	   << "   " << fBuffer[3]
      //      	   << endl;

      detNum[i]  = iBuffer;
      Edeposited = fBuffer[0];

      ElabAB += Edeposited;               // Addback
      Edep[detNum[i]] += Edeposited;      // No addback

    }

    Int_t nDets = 0;

    if(ElabAB > 0) {
      ElabAB = 0; // (We'll recalculate it with simulated thresholds.)
      for(Int_t i = 0; i < nHits; i++){
	Float_t E = Edep[detNum[i]];
	if( E > 0.0 ){

	    // Simulate thresholds
	    if( gRandom->Rndm() 
		< (1.0 + tanh((E-threshPar0[crystalNum[detNum[i]]])/threshPar1[crystalNum[detNum[i]]]))/2.0 ){

	      sim->Fill(  measuredE(E) );            // No addback
	      crystal[crystalNum[detNum[i]]]->Fill(  measuredE(E) );
	      ElabAB += E;                           // Addback

	    }

	    nDets++;
	    Edep[detNum[i]] = 0.0; // Very important to zero after use.
	}
      }
      sim_addback->Fill( measuredE(ElabAB) );      // Addback

      // Photopeak efficiencies
      if(abs(ElabAB - Egamma)<0.01){
	nPhotopeak++;
	for(Int_t i = iFirst; i < iLast; i++){

	  Int_t j;
	  Double_t branch;
	  if(sourceType == "eu152"){
	    j = 0; 
	    branch = Branch[i];
	  }
	  else if(sourceType == "co56") {
	    j = 1;
	    branch = Branch[i];
	  }
	  else if(sourceType == "photopeaks") {
	    j = 2;
	    branch = 1.0;
	  }
	  else {
	    cout << "Error: Unknown source type " << sourceType << endl;
	  }

	  if (abs(ElabAB - Energy[i])<0.01) {
	    EffAB[i] += 1.0/gammasEmitted/branch*gammasPerDecay[j];
	    if(nDets == 1)
	      Eff[i] += 1.0/gammasEmitted/branch*gammasPerDecay[j];
	  }
	}
      }
    }

    nGamma++;
    if(nGamma%100000 == 0) {
      cout << "\r... sorted " << nGamma << " gamma rays ...";
      flush(cout);
    }

  }

  fp.close();

  // Make ring spectra by summing crystal spectra.
  // ring_49: 2, 6, 10, 14

  ring[0]->Sumw2();
  ring[0]->Add(crystal[2], crystal[6],  1.0, 1.0);
  ring[0]->Add(crystal[10], 1.0);
  ring[0]->Add(crystal[14], 1.0);

  // ring_55: 1, 5, 9, 13
  ring[1]->Sumw2();
  ring[1]->Add(crystal[1], crystal[5],  1.0, 1.0);
  ring[1]->Add(crystal[9],  1.0);
  ring[1]->Add(crystal[13], 1.0);

  // ring_64: 3, 7, 11, 15
  ring[2]->Sumw2();
  ring[2]->Add(crystal[3], crystal[7],  1.0, 1.0);
  ring[2]->Add(crystal[11], 1.0);
  ring[2]->Add(crystal[15], 1.0);

  // ring_67: 0, 4, 8, 12
  ring[3]->Sumw2();
  ring[3]->Add(crystal[0], crystal[4], 1.0, 1.0);
  ring[3]->Add(crystal[8], 1.0);
  ring[3]->Add(crystal[12], 1.0);

  // ring_74: 19, 23
  ring[4]->Sumw2();
  ring[4]->Add(crystal[19], crystal[23], 1.0, 1.0);

  // ring_87: 16, 20
  ring[5]->Sumw2();
  ring[5]->Add(crystal[16], crystal[20], 1.0, 1.0);

  // ring_93: 18, 22
  ring[6]->Sumw2();
  ring[6]->Add(crystal[18], crystal[22], 1.0, 1.0);

  // ring_105: 17, 21
  ring[7]->Sumw2();
  ring[7]->Add(crystal[17], crystal[21], 1.0, 1.0);

  cout << "\r... sorted " 
       << nGamma     << " gamma events, and " 
       << nPhotopeak << " photopeak events ...\n" << endl;

  return;

};

void UCGretinaSourceSort(TString simFile){

  cout << "\n... loading Doppler reconstruction and resolution parameters ...\n" << endl;
  loadGretina("Gretina.dat");

  TString rootFile = simFile.Copy();
  rootFile.Replace(rootFile.Index(".out",4,0,0),4,".root",5);
  TFile *rF = new TFile(rootFile,"recreate");

  TString sourceType = simFile.Copy();
  sourceType.Replace(sourceType.Index(".out",4,0,0),4,"",0);

  cout << "... sorting events in " << simFile << " ...\n" << endl;

  loadSim(simFile, sourceType);

  cout << "... writing Spectra to " << rootFile << " ...\n" << endl;

  rF->Write();

  if(sourceType.Contains("eu152") || sourceType.Contains("co56") ||
     sourceType.Contains("photopeaks")){
    TString effFile = simFile.Copy();
    effFile.Replace(effFile.Index(".out",4,0,0),4,"_eff.dat",8);
    ofstream eF;
    eF.open(effFile, ios::trunc);

    cout << "... writing photopeak efficiencies to " << effFile 
	 << " ...\n" << endl;

    for(Int_t i = iFirst;i < iLast;i++)
      eF << std::setw(10) << Energy[i] << "\t" 
	 << Eff[i] << "\t"
	 << EffAB[i] << endl;
  }

  return;

}
