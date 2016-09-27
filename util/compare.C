#include "Riostream.h"
#include "iostream"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFile.h"

GH1D *spectrum1;
GH1D *spectrum2;

void compare(TString file1, TString name1, TString file2, TString name2) {
  compare(file1, name1, file2, name2, 1.0);
}

void compare(TString file1, TString name1, TString file2, TString name2,
	     Double_t factor) {

  TFile *sF1 = new TFile(file1);
  sF1->GetObject(name1, spectrum1);

  TFile *sF2 = new TFile(file2);
  sF2->GetObject(name2, spectrum2); 

  spectrum2->Scale(factor);
  
  spectrum1->SetLineColor(1); // 1 = Black
  spectrum1->GetXaxis()->SetRangeUser(0.,2000.);
  spectrum2->SetLineColor(4); // 4 = Blue	

  TH1F *diff = (TH1F*)spectrum1->Clone("diff");

  diff->Add(spectrum2,-1.0); 

  TCanvas *Spectra = new TCanvas("Spectra", "Spectra", 600, 600);
  Spectra->SetCrosshair();      //Include a crosshair cursor...
  Spectra->ToggleEventStatus(); //...and show its coordinates
  Spectra->Divide(1,2);

  Spectra->cd(1);
  Spectra_1->SetLogy();

  spectrum1->Draw();
  spectrum2->Draw("SAME");

  Spectra->cd(2);
  diff->Draw();

  cout << "Total difference = " << diff->Integral() << endl;
  cout << "                   "
       << diff->Integral()/spectrum1->Integral()*100.
       << " %" << endl;
  
  cout << "Difference > 60 keV = " << diff->Integral(60, 2000) << endl;
  cout << "                   "
       << diff->Integral(60,2000)/spectrum1->Integral()*100.
       << " %" << endl;

  cout << "Difference > 1415 keV = " << diff->Integral(1415, 2000) << endl;
  cout << "                   "
       << diff->Integral(1415, 2000)/spectrum1->Integral()*100.
       << " %" << endl;

}
