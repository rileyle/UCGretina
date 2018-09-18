#include "Riostream.h"
#include "iostream"
#include "TRandom.h"

void makeSpec_GrR() {

  // Rean in simulated gamma-ray spectra========================================

  TH1F* sim[3];
  TString simFileName[3];
  simFileName[0] = "s44_1329_histos.root";
  simFileName[1] = "s44_2150_histos.root";
  simFileName[2] = "s44_2457_histos.root";
  for(Int_t i = 0; i < 3; i++){
    TFile *sF = new TFile(simFileName[i]);
    sF->GetObject("egamdc",sim[i]);
  }
  
  // Make the background =======================================================

  TH1F *bg1 = new TH1F("bg1","",8000,0.,8000.);
  TF1 *exp1 = new TF1("exp1","exp(10.0-1e-2*x)",0.,8000.);
  bg1->FillRandom("exp1",5000);

  TH1F *bg2 = new TH1F("bg2","",8000,0.,8000.);
  TF1 *exp2 = new TF1("exp2","exp(5.0-1e-3*x)",0.,8000.);
  bg2->FillRandom("exp2",1000);

  // Simulate thresholds
  Double_t e = 80.0; 
  Double_t de = 5.0;
  for(Int_t i = 0; i<200; i++){
    //0.5*( 1 + tanh( (Energy - E) / dE ) )
    Double_t thresh = 0.5*( 1.0 + tanh( ((Double_t)i - e) / de ) );
    bg1->SetBinContent(i,(Int_t)(thresh*bg1->GetBinContent(i)));
    bg2->SetBinContent(i,(Int_t)(thresh*bg2->GetBinContent(i)));
  }

  TH1F* spectrum = bg1->Clone("s44spectrum");
  spectrum->Add(bg2,1.0);
  for(Int_t i = 0; i < 3; i++)
    spectrum->Add(sim[i],1.0);

  spectrum->SetTitle("Simulated ^{44}S Measurement");
  
  TFile *rF = new TFile("s44_spectrum.root","recreate");
  spectrum->Write();
  rF->Close();

  return;

}
