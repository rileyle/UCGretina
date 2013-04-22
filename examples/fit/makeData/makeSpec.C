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
  TH1F *bg1 = new TH1F("bg1","",8000,0.,8000.);
  TF1 *exp1 = new TF1("exp1","exp(6.5-6.5e-3*x)",0.,8000.);
  bg1->FillRandom("exp1",1000);

  TH1F *bg2 = new TH1F("bg2","",8000,0.,8000.);
  TF1 *exp2 = new TF1("exp2","exp(5.0-8.5e-4*x)",0.,8000.);
  bg2->FillRandom("exp2",10000);

  // Simulate thresholds
  Double_t e = 80.0; 
  Double_t de = 5.0;
  for(Int_t i = 0; i<200; i++){
    //0.5*( 1 + tanh( (Energy - E) / dE ) )
    Double_t thresh = 0.5*( 1.0 + tanh( ((Double_t)i - e) / de ) );
    bg1->SetBinContent(i,(Int_t)(thresh*bg1->GetBinContent(i)));
    bg2->SetBinContent(i,(Int_t)(thresh*bg2->GetBinContent(i)));
  }

  // Load Geant4 simulations ===================================================
  cout << "\n... loading simulated spectra ...\n" << endl;

  const Int_t nSim = 4;       // Number of simulated gamma-rays
  TString simFileName[nSim];
  simFileName[0] = "example_500_sim_histos.root";
  simFileName[1] = "example_1000_sim_histos.root";
  simFileName[2] = "example_1500_sim_histos.root";
  simFileName[3] = "example_2000_sim_histos.root";

  for(int i=0;i<nSim;i++) {
    TFile *sF = new TFile(simFileName[i]);
    sF->GetObject("egamdc",sim[i]);
  }

  for(int i=1;i<nSim;i++) {
    sim[0]->Add(sim[i],1.0);
  }
  
  sim[0]->Add(bg1,1.0);
  sim[0]->Add(bg2,1.0);

  sim[0]->Draw();

  sim[0]->SetName("example_spectrum");
  sim[0]->SetTitle("example_spectrum");

  TFile *rF = new TFile("example_spectrum.root","recreate");
  sim[0]->Write();
  rF->Close();

  return;

}
