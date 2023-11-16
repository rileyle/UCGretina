void angDist(){

  GH1D *sim[2];
  GH1D *theta[2];

  TString simFileName[2];
  simFileName[0] = "s44_1329_histos.root";
  simFileName[1] = "s44_1329_iso_histos.root";
  
  for(int i=0;i<2;i++) {
      TFile *sF = new TFile(simFileName[i]);
      sF->GetObject("sim/emitted_theta", sim[i]);
      sim[i]->SetDirectory(0);
      sF->Close();
  }

  for(int i=0;i<2;i++) {
      TFile *sF = new TFile(simFileName[i]);
      sF->GetObject("position/theta", theta[i]);
      theta[i]->SetDirectory(0);
      sF->Close();
  }

  TCanvas* c = new TCanvas("", "", 600, 800);
  c->Divide(1, 2);
  
  c->cd(1);
  sim[0]->SetLineColor(kRed);
  sim[0]->Draw();
  sim[0]->SetStats(0);
  sim[0]->SetTitle("");
  sim[0]->GetXaxis()->SetTitle("Emitted #theta (deg)");
  sim[0]->GetYaxis()->SetTitle("Counts/bin");
  sim[0]->GetXaxis()->SetTitleSize(.05);
  sim[0]->GetXaxis()->SetTitleOffset(0.8);
  sim[0]->GetYaxis()->SetTitleSize(.05);
  sim[0]->GetYaxis()->SetTitleOffset(1.0);


  sim[1]->SetLineColor(kBlue);
  sim[1]->Draw("SAME");

  c->cd(2);
  theta[0]->SetLineColor(kRed);
  theta[0]->Draw();
  theta[0]->SetStats(0);
  theta[0]->SetTitle("");
  theta[0]->GetXaxis()->SetTitle("#theta (deg)");
  theta[0]->GetYaxis()->SetTitle("Counts/bin");
  theta[0]->GetXaxis()->SetTitleSize(.05);
  theta[0]->GetXaxis()->SetTitleOffset(0.8);
  theta[0]->GetYaxis()->SetTitleSize(.05);
  theta[0]->GetYaxis()->SetTitleOffset(1.0);
  theta[1]->SetLineColor(kBlue);
  theta[1]->Draw("SAME");
  
}

