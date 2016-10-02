#include "Riostream.h"
#include "iostream"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFile.h"

// These are declared in global scope so that they are accessible to the 
// fit function fitf().
GH1D *sim;       
GH1D *spectrum;
GH1D *backGround;

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

void eu152Fit() {
  eu152Fit("eu152_gammas_histos.root", 4);
}

void eu152Fit(TString simFileName) {
  eu152Fit(simFileName, 4);
}

void eu152Fit(TString simFileName, Int_t rebin) {

  // Load measured spectra =====================================================
  cout << "\nLoading measured spectra ...\n" << endl;

  TString specFileName = "data/run020_cal_histos.root";
  TFile *spF = new TFile(specFileName);
  TH1F* spec;
  spF->GetObject("egam", spec);
  spectrum = (GH1D*)spec->Clone();
  spectrum->Rebin(rebin);
  
  TString bgFileName   = "data/run007_cal_histos.root";
  TFile *bgF = new TFile(bgFileName);
  TH1F* back;
  bgF->GetObject("egam", back);
  backGround = (GH1D*)back->Clone();
  backGround->Rebin(rebin);
  
  // Load Geant4 simulations ===================================================
  cout << "\n... loading simulated spectra ...\n" << endl;

  const Int_t nSim = 1;       // Number of simulated gamma-rays
  
  //  TString simFileName = "eu152_gammas_histos.root";
  TFile *sF = new TFile(simFileName);
  sF->GetObject("energy/energy_gaus", sim);
  sim->Rebin(rebin);
  
  // Check binning.
  cout << "\n   Measured spectrum has " << spectrum->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   Simulations have "      << sim->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   These should match for meaninful absolute values of fit parameters.\n" << endl;


  // Fit =======================================================================
  cout << "\n... fitting ...\n" << endl;

  const Int_t nPar = nSim+3;  // Number of fit parameters
  Float_t bgMinE  = 1800.;
  Float_t bgMaxE = 3000.;

  Float_t fitMinE = 10.;
  Float_t fitMaxE = 3000.;

  // Initialize fit parameter array
  Double_t par[nPar];
  par[0]  = 0.1;      // Measured random background
  par[1]  = 10.0;     // Simulation
  par[2]  = 85.;      // Threshold parameter 1
  par[3]  =  5.;      // Threshold parameter 2

  TF1 *f1 = new TF1("f1", fitf, bgMinE, bgMaxE, nPar);
  f1->SetLineColor(kBlue);
  f1->SetNpx(10000);        // Number of points used to draw the function

  f1->SetParameters(par);

  // Fit the room background component.
  f1->FixParameter(1, 10.);
  f1->FixParameter(2, 85.);
  f1->FixParameter(3, 5.);

  spectrum->GetXaxis()->SetRangeUser(0.,3000.);
  GH1D *diff = (GH1D*)spectrum->Clone("diff");
  GH1D * sqrtup = (GH1D*)spectrum->Clone("Error+");
  GH1D * sqrtdn = (GH1D*)spectrum->Clone("Error-");

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 1000000 );
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations( 100000 );

  spectrum->Fit("f1","LR");
  
  f1->FixParameter(0, f1->GetParameter(0));

  // Fit the full spectrum.
  f1->SetRange(fitMinE, fitMaxE);

  f1->ReleaseParameter(1);
  f1->ReleaseParameter(2);
  f1->ReleaseParameter(3);
  spectrum->Fit("f1","LR");

  // Discrepancy ===============================================================
  
  // Create the difference spectrum
  for(int i=0;i<diff->GetNbinsX();i++) {

    Float_t E = diff->GetBinCenter(i);

    if(E>fitMinE && E<fitMaxE) {
      diff->SetBinContent(i, f1->Eval(E) - diff->GetBinContent(i));
      sqrtup->SetBinContent(i,  sqrt(spectrum->GetBinContent(i)
				     +f1->Eval(E)));
      sqrtdn->SetBinContent(i, -sqrt(spectrum->GetBinContent(i)
				     +f1->Eval(E)));

    } else {
      diff->SetBinContent(i, 0.);
      sqrtup->SetBinContent(i, 0.);
      sqrtdn->SetBinContent(i, 0.);
    }

  }

  // Display ===================================================================
  
  TCanvas *Spectra = new TCanvas("Spectra", "Spectra", 800, 800);
  Spectra->SetCrosshair();      //Include a crosshair cursor...
  Spectra->ToggleEventStatus(); //...and show its coordinates
  Spectra->Divide(1,2);

  Spectra->cd(1);
  Spectra_1->SetLogy();

  spectrum->SetStats(kFALSE);
  spectrum->SetLineColor(kBlack);
  spectrum->GetXaxis()->SetTitle("Energy (keV)");
  spectrum->GetYaxis()->SetTitle(Form("Counts/(%d keV)", rebin));
  spectrum->Draw();

  Spectra->cd(2);

  diff->SetTitle("Discrepancy");
  diff->SetStats(kFALSE);
  diff->GetXaxis()->SetTitle("Energy (keV)");
  diff->GetYaxis()->SetTitle(Form("Counts/(%d keV)", rebin));
  diff->SetFillColor(kBlue);
  diff->SetLineColor(kBlue);
  diff->GetYaxis()->SetRangeUser(-5000.,5000.);
  diff->Draw();

  sqrtup->SetFillColor(kGray);
  sqrtup->SetLineColor(kGray);
  sqrtup->Draw("SAME");
  sqrtdn->SetFillColor(kGray);
  sqrtdn->SetLineColor(kGray);
  sqrtdn->Draw("SAME");
  diff->Draw("SAME");

  cout << "Total difference = " << diff->Integral() << endl;

}
