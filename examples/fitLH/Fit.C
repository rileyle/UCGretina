#include "Riostream.h"
#include "iostream"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFile.h"

// These are declared in global scope so that they are accessible to the 
// fit function fitf().
TH1F *sim[12];       
TH1F *spectrum;
TH1F *backGround;

// Fit function
Double_t fitf(Double_t *v,Double_t *par) {

  // Get bin contents at v[0]
  Double_t bg = 0; // Remove and uncomment next line for random background
  //  Double_t bg    = backGround->GetBinContent(backGround->FindBin(v[0]));
  Double_t sim0  = sim[0]->GetBinContent(sim[0]->FindBin(v[0]));
  Double_t sim1  = sim[1]->GetBinContent(sim[1]->FindBin(v[0]));
  Double_t sim2  = sim[2]->GetBinContent(sim[2]->FindBin(v[0]));
  Double_t sim3  = sim[3]->GetBinContent(sim[3]->FindBin(v[0]));

  // Return a linear combination of two exponentials, the measured 
  // background and the simulated gamma rays.
  Double_t value = exp(par[0]+par[1]*v[0]) + exp(par[2]+par[3]*v[0])
    +  par[4]*bg
    +  par[5]*sim0 +  par[6]*sim1 +  par[7]*sim2  +  par[8]*sim3;

  return value;

}

void fitSpectrum() {

  // Load measured spectra =====================================================
  cout << "\nLoading measured spectrum ...\n" << endl;

  TString specFileName = "example_spectrum.root";
  TFile *spF = new TFile(specFileName);
  spF->GetObject("example_spectrum", spectrum);
  spectrum->Rebin(8);

  // To include random background, need a room background spectrum Doppler
  // corrected with the same beta is the measured spectrum.

  //  TString bgFileName   = "??.root";  
  //  TFile *bgF = new TFile(bgFileName);
  //  bgF->GetObject("egamdc", backGround);

  // Load Geant4 simulations ===================================================
  cout << "\n... loading simulated spectra ...\n" << endl;

  const Int_t nSim = 4;       // Number of simulated gamma-rays
  TString simFileName[nSim];
  simFileName[0] = "example_500";
  simFileName[1] = "example_1000";
  simFileName[2] = "example_1500";
  simFileName[3] = "example_2000";

  for(int i=0;i<nSim;i++) {
    simFileName[i] += "_sim_histos.root";
    TFile *sF = new TFile(simFileName[i]);
    sF->GetObject("egamdc",sim[i]);
    sim[i]->Rebin(8);
  }

  // Check binning.
  cout << "\n   Measured spectrum has " << spectrum->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   Simulations have "      << sim[0]->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   These should match for meaninful absolute values of fit parameters.\n" << endl;

  // Fit =======================================================================
  cout << "\n... fitting ...\n" << endl;

  const Int_t nPar = nSim+5;  // Number of fit parameters
  Float_t fitMinE = 110.;
  Float_t fitMaxE = 4000.;

  // Initialize fit parameter array
  Double_t par[1000];
  par[0]  = 1;                // Exponential 1
  par[1]  = -1E-03;
  par[2]  = 1;                // Exponential 2
  par[3]  = -1E-04;
  par[4]  = 0;                // Measured random background
  par[5]  = 0.001;             //  500 keV
  par[6]  = 0.002;             // 1000 keV
  par[7]  = 0.005;             // 1500 keV
  par[8]  = 0.01;              // 2000 keV

  TF1 *f1 = new TF1("f1",fitf,fitMinE,fitMaxE,nPar);
  f1->SetLineColor(4); // 4 = Blue
  f1->SetNpx(1000);     // Number of points used to draw the function
  f1->SetLineWidth(2);

  f1->SetParameters(par);
  f1->FixParameter(4,0); // No actual measured background here ...

  spectrum->GetXaxis()->SetRangeUser(0.,2500.);
  TH1F *diff = (TH1F*)spectrum->Clone("diff");

  spectrum->Fit("f1","LR");

  // Create the difference spectrum
  for(int i=0;i<diff->GetNbinsX();i++) {

    Float_t E = diff->GetBinCenter(i);

    if(E>fitMinE && E<fitMaxE) {
      diff->SetBinContent(i, diff->GetBinContent(i) - f1->Eval(E));
    } else {
      diff->SetBinContent(i, 0.);
    }

  }

  TCanvas *Spectra = new TCanvas("Spectra", "Spectra", 600, 600);
  Spectra->SetCrosshair();      //Include a crosshair cursor...
  Spectra->ToggleEventStatus(); //...and show its coordinates
  Spectra->Divide(1,2);

  Spectra->cd(1);
  Spectra_1->SetLogy();

  spectrum->Draw();

  Spectra->cd(2);
  //diff->GetYaxis()->SetRangeUser(-100.,100.);
  diff->Draw();

}
