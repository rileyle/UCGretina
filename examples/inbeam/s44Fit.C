
// These must be declared in global scope so that they are accessible to the 
// fit function fitf().
GH1D *sim[12];       
GH1D *spectrum;
GH1D *backGround;

// Fit function
Double_t fitf(Double_t *v,Double_t *par) {

  // Get bin contents at v[0]
  Double_t bg = 0; // Remove and uncomment next line for random background
  //  Double_t bg    = backGround->GetBinContent(backGround->FindBin(v[0]));
  Double_t sim0  = sim[0]->GetBinContent(sim[0]->FindBin(v[0]));
  Double_t sim1  = sim[1]->GetBinContent(sim[1]->FindBin(v[0]));
  Double_t sim2  = sim[2]->GetBinContent(sim[2]->FindBin(v[0]));

  // Return a linear combination of two exponentials, the measured 
  // background and the simulated gamma rays.
  Double_t value = exp(par[0]+par[1]*v[0]) + exp(par[2]+par[3]*v[0])
    +  par[4]*bg
    +  par[5]*sim0 +  par[6]*sim1 +  par[7]*sim2;

  return value;

}

void s44Fit() {
  s44Fit(8);
}

void s44Fit(Int_t rebin) {

  // Load measured spectra =====================================================
  cout << "\nLoading measured spectrum ...\n" << endl;

  TString specFileName = "s44_spectrum.root";
  TFile *spF = new TFile(specFileName);
  spF->GetObject("s44spectrum", spectrum);
  spectrum->Rebin(rebin);

  // To include random background, need a room background spectrum Doppler
  // corrected with the same beta is the measured spectrum.

  //  TString bgFileName   = "??.root";  
  //  TFile *bgF = new TFile(bgFileName);
  //  bgF->GetObject("egamdc", backGround);

  // Load Geant4 simulations ===================================================
  cout << "\n... loading simulated spectra ...\n" << endl;

  const Int_t nSim = 3;       // Number of simulated excited states
  TString simFileName[nSim];
  simFileName[0] = "s44_1329_histos.root";
  simFileName[1] = "s44_2150_histos.root";
  simFileName[2] = "s44_2457_histos.root";

  for(int i=0;i<nSim;i++) {
    TFile *sF = new TFile(simFileName[i]);
    sF->GetObject("energy/dop_4169_gaus",sim[i]);
    sim[i]->Rebin(rebin);
  }

  // Check binning.
  cout << "\n   Measured spectrum has " << spectrum->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   Simulations have "      << sim[0]->GetBinWidth(1) 
       << " keV/channel" << endl;
  cout << "   These must match for meaninful absolute values of fit parameters.\n" << endl;

  // Fit =======================================================================
  cout << "\n... fitting ...\n" << endl;

  const Int_t nPar = nSim+5;  // Number of fit parameters
  Double_t fitMinE = 110.;
  Double_t fitMaxE = 3000.;

  // Initialize fit parameter array
  Double_t par[1000];
  par[0]  = 1.;              // Exponential 1
  par[1]  = -1E-02;
  par[2]  = 1.;              // Exponential 2
  par[3]  = -1E-03;
  par[4]  = 0;                // Measured random background
  par[5]  = 0.0001;           // 1329 keV
  par[6]  = 0.0001;           // 2150 keV
  par[7]  = 0.0001;           // 2457 keV

  TF1 *f1 = new TF1("f1",fitf,fitMinE,fitMaxE,nPar);
  f1->SetLineColor(kBlue);
  f1->SetNpx(1000);     // Number of points used to draw the function
  f1->SetLineWidth(2);

  f1->SetParameters(par);
  f1->FixParameter(4,0); // No actual measured background here ...

  f1->SetParLimits(0, 0, 100);
  f1->SetParLimits(1, -0.1, -0.005);
  f1->SetParLimits(2, 0, 100);
  f1->SetParLimits(3, -0.005, 0);
  
  spectrum->GetXaxis()->SetRangeUser(0., fitMaxE);
  GH1D *diff    = (GH1D*)spectrum->Clone("diff");
  GH1D * sqrtup = (GH1D*)spectrum->Clone("Error+");
  GH1D * sqrtdn = (GH1D*)spectrum->Clone("Error-");

  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Migrad");
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 1000000 );
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations( 100000 );
   
  spectrum->Fit("f1","LR");

  // Discrepancy ===============================================================
  
  // Create the difference spectrum and channel error spectra
  for(int i=0;i<diff->GetNbinsX();i++) {

    Double_t E = diff->GetBinCenter(i);

    if(E>fitMinE && E<fitMaxE) {
      diff->SetBinContent(i, diff->GetBinContent(i) - f1->Eval(E));
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
  
  TCanvas *Spectra = new TCanvas("Spectra", "Spectra", 600, 600);
  Spectra->SetCrosshair();      //Include a crosshair cursor...
  Spectra->ToggleEventStatus(); //...and show its coordinates
  Spectra->Divide(1,2);

  Spectra->cd(1);
  Spectra_1->SetLogy();

  spectrum->SetLineColor(kBlack);
  spectrum->Draw();

  Spectra->cd(2);

  diff->SetTitle("Discrepancy");
  diff->SetFillColor(kBlue);
  diff->SetLineColor(kBlue);

  diff->Draw();
  sqrtup->SetFillColor(kGray);
  sqrtup->SetLineColor(kGray);
  sqrtup->Draw("SAME");
  sqrtdn->SetFillColor(kGray);
  sqrtdn->SetLineColor(kGray);
  sqrtdn->Draw("SAME");
  diff->Draw("SAME");

  cout << "\nPopulation of excited states of 44S:" << endl;
  cout << "   1329 keV: " << f1->GetParameter(5)*1e6 << " +/- "
       << f1->GetParError(5)*1e6 << endl;
  cout << "   2150 keV: " << f1->GetParameter(6)*1e6 << " +/- "
       << f1->GetParError(6)*1e6 << endl;
  cout << "   2457 keV: " << f1->GetParameter(7)*1e6 << " +/- "
       << f1->GetParError(7)*1e6 << endl;

}
