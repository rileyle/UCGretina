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

// Doppler-corrected energy with simulated tracking resolution
Float_t DopplerE(Float_t E, Float_t x, Float_t y, Float_t z){

  x += gRandom->Gaus(0,2.0);
  y += gRandom->Gaus(0,2.0);
  z += gRandom->Gaus(0,2.0);
  z -= zoffset;

  Float_t cosTheta = z/sqrt(x*x + y*y + z*z);
	
  Float_t gamma = 1/sqrt(1-beta*beta);

  return E*(1-beta*cosTheta)*gamma; 
}

// Read simulation data and sort into a Doppler-corrected histogram
void loadSim(TString fileName) {

  char line[1000];
  Float_t fBuffer[4];
  Int_t   iBuffer;
  TString sBuffer;

  Float_t Edep[1000] = {0.0};
  Float_t xHit[1000] = {0.0};
  Float_t yHit[1000] = {0.0};
  Float_t zHit[1000] = {0.0};
  Int_t   detNum[200] = {0};

  ifstream fp;
  fp.open(fileName);

  Int_t nChannels     = 1024;
  Int_t keVperChannel = 8;    // Match binning of the spectrum to fit
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
      x          = fBuffer[1];
      y          = fBuffer[2];
      z          = fBuffer[3];

      // Addback
      ElabAB += Edeposited;
      //    Set the hit position = position of largest energy deposit.
      if(Edeposited > EdepMax){
	EdepMax = Edeposited;
	xs = x;
	ys = y;
	zs = z;
      }

      // No addback
      //   Hit position = first hit in crystal
      if(!Edep[detNum[i]] > 0.){
	xHit[detNum[i]] = x;
	yHit[detNum[i]] = y;
	zHit[detNum[i]] = z;
      }
      Edep[detNum[i]] += Edeposited;

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

	    // No addback
	    sim->Fill( DopplerE( measuredE(E),
				 xHit[detNum[i]],
				 yHit[detNum[i]],
				 zHit[detNum[i]]) );
	    // No addback
	    crystal[crystalNum[detNum[i]]]->Fill( DopplerE( measuredE(E),     
							    xHit[detNum[i]],
							    yHit[detNum[i]],
							    zHit[detNum[i]]) );

	    // Addback			
	    ElabAB += E;

	  }

	  nDets++;
	  Edep[detNum[i]] = 0.0; // Very important to zero after use.
	  xHit[detNum[i]] = 0.;
	  yHit[detNum[i]] = 0.;
	  zHit[detNum[i]] = 0.;

	}
      }

      // Addback
      sim_addback->Fill( DopplerE(measuredE(ElabAB), xs, ys, zs) );

      // Photopeak events
      if( abs(Egamma - ElabAB) < 0.02 )
	nPhotopeak++;
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

void UCGretinaSort(TString simFile){

  cout << "\n... loading Doppler reconstruction and resolution parameters ...\n" << endl;
  loadGretina("Gretina.dat");

  TString rootFile = simFile.Copy();
  rootFile.Replace(rootFile.Index(".out",4,0,0),4,".root",5);
  TFile *rF = new TFile(rootFile,"recreate");

  TString sourceType = simFile.Copy();
  sourceType.Replace(sourceType.Index(".out",4,0,0),4,"",0);

  cout << "... sorting events in " << simFile << " ...\n" << endl;

  loadSim(simFile);

  cout << "... writing Spectra to " << rootFile << " ...\n" << endl;

  rF->Write();

  return;

}
