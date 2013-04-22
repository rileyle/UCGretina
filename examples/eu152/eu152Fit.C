#include "Riostream.h"
#include "iostream"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFile.h"

// These are declared in global scope so that they are accessible to the 
// fit function fitf().
TH1F *sim;       
TH1F *spectrum;
TH1F *backGround;

// Fit function
Double_t fitf(Double_t *v,Double_t *par) {

  // Get bin contents at v[0]
  Double_t bg    = backGround->GetBinContent(backGround->FindBin(v[0]));
  Double_t sim0  = sim->GetBinContent(sim->FindBin(v[0]));

  // Threshold
  sim0 = sim0*(1.0 + tanh((v[0]-par[2])/par[3]))/2.0;

  // Return a linear combination of the measured 
  // background and the simulated gamma rays.
  Double_t value = par[0]*bg + par[1]*sim0;

  return value;

}

void fitSpectrum() {
  fitSpectrum("egam");
}

void fitSpectrum(TString histoName) {

  // Load measured spectra =====================================================
  cout << "\nLoading measured spectra ...\n" << endl;

  TString specFileName = "data/run020_cal_histos.root";
  TFile *spF = new TFile(specFileName);
  spF->GetObject(histoName,spectrum);

  TString bgFileName   = "data/run007_cal_histos.root";
  TFile *bgF = new TFile(bgFileName);
  bgF->GetObject(histoName,backGround);

  // Load Geant4 simulations ===================================================
  cout << "\n... loading simulated spectra ...\n" << endl;

  const Int_t nSim = 1;       // Number of simulated gamma-rays
  
  TString simFileName = "eu152_sim_histos.root";
  TFile *sF = new TFile(simFileName);
  sF->GetObject(histoName,sim);

  // Check binning.
  cout << "\n   Measured spectrum has " << spectrum->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   Simulations have "      << sim->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   These should match for meaninful absolute values of fit parameters.\n" << endl;


  // Fit =======================================================================
  cout << "\n... fitting ...\n" << endl;

  const Int_t nPar = nSim+3;  // Number of fit parameters
  Float_t bgMinE  = 1540.;
  Float_t bgMaxE = 2500.;

  Float_t fitMinE = 10.;
  Float_t fitMaxE = 4000.;

  // Initialize fit parameter array
  Double_t par[1000];
  par[0]  = 0.1;      // Measured random background
  par[1]  = 2.0;             // Simulation
  par[2]  = 50.;           // Threshold parameter 1
  par[3]  =  5.;           // Threshold parameter 2

  TF1 *f1 = new TF1("f1",fitf,bgMinE,bgMaxE,nPar);
  f1->SetLineColor(4); // 4 = Blue
  f1->SetNpx(1000);    // Number of points used to draw the function

  f1->SetParameters(par);
  //  f1->FixParameter(0, 0.109163); // Run 081 (Based on 1540-4000 keV fit)
  f1->FixParameter(2,0.0);
  f1->FixParameter(3,0.001);

  spectrum->GetXaxis()->SetRangeUser(0.,4000.);
  TH1F *diff = (TH1F*)spectrum->Clone("diff");

  // Set the background.
  spectrum->Fit("f1","LRME");
  
  f1->FixParameter(0,f1->GetParameter(0));

  // Fit the full spectrum.
  f1->SetRange(fitMinE, fitMaxE);
  spectrum->Fit("f1","LRME");
  
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

  spectrum->SetLineColor(1);
  spectrum->Draw();

  Spectra->cd(2);
  diff->SetLineColor(1);
  diff->GetYaxis()->SetRangeUser(-2000.,2000.);
  diff->Draw();

  cout << "Total difference = " << diff->Integral() << endl;

}
