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

// Load a spectrum in NSCL ascii format into a root TH1F histogram,
// and return a pointer to the histogram.
TH1F *loadSpec(TString fileName) {

  char line[1000];

  ifstream fp;
  fp.open(fileName);

  char name[30];
  int nChannels;
  int lo,hi;

  //Header line 1 : "spectrum name"  (nChannels)
  fp.getline(line,1000);
  sscanf(line,"%s (%d)",&name, &nChannels);

  //Header line 2 : Date string
  fp.getline(line,1000);

  //Header line 3: the number 2!?
  fp.getline(line,1000);

  //Header line 4: dimension type
  fp.getline(line,1000);

  //Header line 5: ("Title")
  fp.getline(line,1000);

  //Header line 6: (lo hi)
  fp.getline(line,1000);
  sscanf(line,"(%d %d)",&lo, &hi);

  //Header line 7: ------------------------------------
  fp.getline(line,1000);

  cout << name << " " << nChannels << " " << lo << " " << hi << endl;

  TH1F *data = new TH1F("data","", nChannels, float(lo), float(hi+1));

  int channel;
  float content;
  while(1)
    {
      if(!(fp.good())) break;

      fp.getline(line,1000);
      sscanf(line,"(%d) %f",&channel,&content);
      data->SetBinContent(channel,content);

    }

  fp.close();

  return data;
}

// Fit function
Double_t fitf(Double_t *v,Double_t *par) {

  // Calculate calibrated energy (for simulations, not background!)
  Double_t eCal = par[2] + par[3]*v[0];

  // Get bin contents at v[0]
  Double_t bg    = backGround->GetBinContent(backGround->FindBin(v[0]));
  Double_t sim0  = sim->GetBinContent(sim->FindBin(eCal));

  // Threshold
  sim0 = sim0*(1.0 + tanh((v[0]-par[4])/par[5]))/2.0;

  // Return a linear combination of the measured 
  // background and the simulated gamma rays.
  Double_t value = par[0]*bg + par[1]*sim0; // + exp(par[4]+par[5]*v[0]);

  return value;

}

void fitSpectrum(TString simName) {

  // Load measured spectra =====================================================
  cout << "\nLoading measured spectra ...\n" << endl;

  TString specFile = "eu152_run087_cal.asc";
  TString bgFile   = "bg_run081_cal.asc";

  spectrum   = loadSpec(specFile);
  //spectrum->Rebin(4);

  backGround = loadSpec(bgFile);
  //backGround->Rebin(4);

  // Load Geant4 simulations ===================================================
  cout << "\n... loading simulated spectra ...\n" << endl;

  const Int_t nSim = 1;       // Number of simulated gamma-rays
  
  TString simFile = simName.Copy();
  simFile += ".root";
  TFile *sF = new TFile(simFile);
  sF->GetObject(simName,sim);

  // Check binning.
  cout << "\n   Measured spectrum has " << spectrum->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   Simulations have "      << sim->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   These should match for meaninful absolute values of fit parameters.\n" << endl;


  // Fit =======================================================================
  cout << "\n... fitting ...\n" << endl;

  const Int_t nPar = nSim+5;  // Number of fit parameters
  Float_t fitMinE = 10.;
  Float_t fitMaxE = 3000.;

  // Initialize fit parameter array
  Double_t par[1000];
  par[0]  = 0.1;      // Measured random background
  par[1]  = 2.0;             // Simulation
  par[2]  = 0;             // Calibration offset
  par[3]  = 1;             // Calibration slope
  par[4]  = 80.;           // Threshold parameter 1
  par[5]  = 20.;           // Threshold parameter 2

  TF1 *f1 = new TF1("f1",fitf,fitMinE,fitMaxE,nPar);
  f1->SetLineColor(4); // 4 = Blue
  f1->SetNpx(1000);    // Number of points used to draw the function

  f1->SetParameters(par);
  //  f1->FixParameter(0, 0.109163); // Run 081 (Based on 1540-4000 keV fit)
  f1->FixParameter(2,0.0);
  f1->FixParameter(3,1.0);
  f1->FixParameter(4,0.0);
  f1->FixParameter(5,0.001);

  spectrum->GetXaxis()->SetRangeUser(0.,3000.);
  TH1F *diff = (TH1F*)spectrum->Clone("diff");

  spectrum->Fit("f1","R");

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
