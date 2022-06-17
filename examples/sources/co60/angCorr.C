void angCorr(){

  GH1D *delta;
  GH1D *isoDelta;
  GH1D *mDelta;
  GH1D *misoDelta;
  GH2D *deltamDelta;
  
  TString fileName    = "co60_histos.root";
  TString isoFileName = "iso_histos.root";

  TFile *sF = new TFile("co60_histos.root");
  sF->GetObject("sim/emitted_delta", delta);
  sF->GetObject("position/delta", mDelta);
  sF->GetObject("position/cosTheta_mCosTheta", deltamDelta);
  
  TFile *iF = new TFile("iso_histos.root");
  iF->GetObject("sim/emitted_delta", isoDelta);
  iF->GetObject("position/delta", misoDelta);

  Int_t rebin = 20;
  if(delta->Integral()>=1e7)
    rebin = 10;
  
  mDelta->Rebin(rebin);
  misoDelta->Rebin(rebin);
  
  TCanvas *c1 = new TCanvas("Angular Correlation", "Angular Correlation",
			    1200, 800);
  c1->SetCrosshair();      //Include a crosshair cursor...
  c1->ToggleEventStatus(); //...and show its coordinates
  c1->Divide(2,2);
  
  c1->cd(1);
  mDelta->SetLineColor(kRed);
  mDelta->SetStats(kFALSE);
  mDelta->SetTitle("");
  mDelta->GetXaxis()->SetTitle("\"Measured\" cos (#theta)");
  mDelta->GetYaxis()->SetTitle("\"Measured\" Counts/bin");
  mDelta->GetXaxis()->SetTitleSize(.05);
  mDelta->GetXaxis()->SetTitleOffset(1.0);
  mDelta->GetYaxis()->SetTitleSize(.05);
  mDelta->GetYaxis()->SetTitleOffset(1.0);
  mDelta->Draw();
  misoDelta->SetLineColor(kBlue);
  misoDelta->Draw("SAME");

  c1->cd(3);
  mDelta->Sumw2();
  misoDelta->Sumw2();
  GH1D *W = (GH1D*)mDelta->Clone();
  W->SetStats(kFALSE);
  W->Divide(misoDelta);
  W->SetMarkerStyle(20);
  W->Draw("PE1");
  W->SetTitle("");
  W->GetXaxis()->SetTitle("\"Measured\" cos (#theta)");
  W->GetYaxis()->SetTitle("W(#theta)");
  //  W->Draw();

  c1->cd(2);
  delta->SetLineColor(kRed);
  delta->SetStats(kFALSE);
  delta->SetTitle("");
  delta->GetXaxis()->SetTitle("Emitted cos (#theta)");
  delta->GetYaxis()->SetTitle("Emitted Counts/bin");
  delta->GetXaxis()->SetTitleSize(.05);
  delta->GetXaxis()->SetTitleOffset(1.0);
  delta->GetYaxis()->SetTitleSize(.05);
  delta->GetYaxis()->SetTitleOffset(1.1);
  delta->Draw();
  isoDelta->SetLineColor(kBlue);
  isoDelta->Draw("SAME");
  
  c1->cd(4);
  deltamDelta->Draw("col");
  deltamDelta->SetStats(kFALSE);
  deltamDelta->SetTitle("");
  deltamDelta->GetXaxis()->SetTitle("Emitted cos (#theta)");
  deltamDelta->GetYaxis()->SetTitle("\"Measured\" cos (#theta)");
  deltamDelta->GetXaxis()->SetTitleSize(.05);
  deltamDelta->GetXaxis()->SetTitleOffset(1.0);
  deltamDelta->GetYaxis()->SetTitleSize(.05);
  deltamDelta->GetYaxis()->SetTitleOffset(1.0);

}

