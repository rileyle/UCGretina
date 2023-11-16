void pol(){

  const Int_t nSim = 2;
  GH1D *sim[nSim];       
  TString simFileName[nSim];
  simFileName[0] = "s44_1329_histos.root";
  simFileName[1] = "s44_1329_nopol_histos.root";

  for(int i=0;i<nSim;i++) {
    TFile *sF = new TFile(simFileName[i]);
    sF->GetObject("polarization/xi", sim[i]);
    sim[i]->Rebin(10);
  }
  
  TCanvas* c = new TCanvas();
  
  sim[0]->SetLineColor(kRed);
  sim[0]->Draw();
  sim[0]->SetStats(0);
  sim[0]->SetTitle("");
  sim[0]->GetXaxis()->SetTitle("#xi (deg)");
  sim[0]->GetYaxis()->SetTitle("Counts/bin");
  sim[0]->GetXaxis()->SetTitleSize(.05);
  sim[0]->GetXaxis()->SetTitleOffset(0.8);
  sim[0]->GetYaxis()->SetTitleSize(.05);
  sim[0]->GetYaxis()->SetTitleOffset(1.0);
  sim[1]->SetLineColor(kBlue);
  sim[1]->Draw("SAME");

}
