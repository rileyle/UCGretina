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

  TH1F *data = new TH1F("data","", nChannels, float(lo), float(hi));

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

  // Get bin contents at v[0]
  Double_t bg    = backGround->GetBinContent(backGround->FindBin(v[0]));
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
  cout << "\nLoading measured spectra ...\n" << endl;

  TString specFile = "example_spectrum.asc";
  TString bgFile   = "bg_dop_all.asc";

  spectrum   = loadSpec(specFile);

  backGround = loadSpec(bgFile);
  backGround->Rebin(4);

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
  par[5]  = 0.01;             //  500 keV
  par[6]  = 0.02;             // 1000 keV
  par[7]  = 0.05;             // 1500 keV
  par[8]  = 0.1;              // 2000 keV

  TF1 *f1 = new TF1("f1",fitf,fitMinE,fitMaxE,nPar);
  f1->SetLineColor(4); // 4 = Blue
  f1->SetNpx(1000);     // Number of points used to draw the function
  f1->SetLineWidth(2);

  f1->SetParameters(par);
  f1->FixParameter(4,0); // No actual measured background here ...

  spectrum->GetXaxis()->SetRangeUser(0.,2500.);
  TH1F *diff = (TH1F*)spectrum->Clone("diff");

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

  spectrum->Draw();

  Spectra->cd(2);
  //diff->GetYaxis()->SetRangeUser(-100.,100.);
  diff->Draw();

}
