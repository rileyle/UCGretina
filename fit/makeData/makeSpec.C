#include "Riostream.h"
#include "iostream"
#include "TVirtualFitter.h"
#include "TRandom.h"

// These are declared in global scope so that they are accessible to the 
// fit function fitf().
TH1F *sim[12];       
TH1F *spectrum;
TH1F *backGround;

void makeSpec() {

  // Make the background =======================================================
  TH1F *bg1 = new TH1F("bg1","",1024,0.,8192.);
  TF1 *exp1 = new TF1("exp1","exp(6.5-6.5e-3*x)",0.,8192.);
  bg1->FillRandom("exp1",1000);

  TH1F *bg2 = new TH1F("bg2","",1024,0.,8192.);
  TF1 *exp2 = new TF1("exp2","exp(5.0-8.5e-4*x)",0.,8192.);
  bg2->FillRandom("exp2",10000);

  // Load Geant4 simulations ===================================================
  cout << "\n... loading simulated spectra ...\n" << endl;

  const Int_t nSim = 4;       // Number of simulated gamma-rays
  TString simName[nSim];
  simName[0] = "example_500";
  simName[1] = "example_1000";
  simName[2] = "example_1500";
  simName[3] = "example_2000";

  for(int i=0;i<nSim;i++) {
    TString simFile = simName[i].Copy();
    simFile += ".root";
    TFile *sF = new TFile(simFile);
    sF->GetObject(simName[i],sim[i]);
  }

  for(int i=1;i<nSim;i++) {
    sim[0]->Add(sim[i],1.0);
  }
  
  sim[0]->Add(bg1,1.0);
  sim[0]->Add(bg2,1.0);

  sim[0]->Draw();

  sim[0]->SetName("Example Spectrum");

  TFile *rF = new TFile("example_spectrum.root","recreate");
  sim[0]->Write();
  rF->Close();

  return;

}
