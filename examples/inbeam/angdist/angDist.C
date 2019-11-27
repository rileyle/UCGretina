void angDist(){

  const Int_t nSim = 2;
  GH1D *sim[nSim];       
  TString simFileName[nSim];
  simFileName[0] = "s44_1329_histos.root";
  simFileName[1] = "s44_1329_iso_histos.root";

  for(int i=0;i<nSim;i++) {
    TFile *sF = new TFile(simFileName[i]);
    sF->GetObject("sim/emitted_theta",sim[i]);
  }

  sim[0]->SetLineColor(kRed);
  sim[0]->Draw();
  sim[1]->SetLineColor(kBlue);
  sim[1]->Draw("SAME");
    
}

