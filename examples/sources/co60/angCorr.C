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
			    800, 800);
  c1->SetCrosshair();      //Include a crosshair cursor...
  c1->ToggleEventStatus(); //...and show its coordinates
  c1->Divide(2,2);

  c1->cd(1);
  mDelta->SetLineColor(kRed);
  mDelta->Draw();
  misoDelta->SetLineColor(kBlue);
  misoDelta->Draw("SAME");

  c1->cd(3);
  mDelta->Sumw2();
  misoDelta->Sumw2();
  GH1D *W = (GH1D*)mDelta->Clone();
  W->Divide(misoDelta);
  W->SetMarkerStyle(20);
  W->Draw("PE1");
  //  W->Draw();

  c1->cd(2);
  delta->SetLineColor(kRed);
  delta->Draw();
  isoDelta->SetLineColor(kBlue);
  isoDelta->Draw("SAME");
  
  c1->cd(4);
  deltamDelta->Draw("col");
}

