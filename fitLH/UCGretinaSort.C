#include "Riostream.h"
#include "iostream"
#include "TRandom.h"
#include "TFile.h"

TH1F *sim;
TH1F *sim_addback;

Float_t beta, zoffset;
Float_t sigmaPar0, sigmaPar1, sigmaPar2, sigmaPar3;
Float_t ecalPar0, ecalPar1;

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
      if(Edep[detNum[i]] > 0.){
	xHit[detNum[i]] = x;
	yHit[detNum[i]] = y;
	zHit[detNum[i]] = z;
      }
      Edep[detNum[i]] += Edeposited;

    }

    Int_t nDets = 0;

    if(ElabAB > 0) {
      // Addback
      sim_addback->Fill( DopplerE(measuredE(ElabAB), xs, ys, zs) );

      // No addback
      for(Int_t i = 0; i < nHits; i++){
	if(Edep[detNum[i]] > 0.0){
	  sim->Fill( DopplerE(measuredE(Edep[detNum[i]]),
			      xHit[detNum[i]],
			      yHit[detNum[i]],
			      zHit[detNum[i]]) );
	  nDets++;
	  Edep[detNum[i]] = 0.0; // Very important to zero after use.
	  xHit[detNum[i]] = 0.;
	  yHit[detNum[i]] = 0.;
	  zHit[detNum[i]] = 0.;
	}
      }

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
