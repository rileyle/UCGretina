void energies(){
  energies(31);
}

void energies(Int_t mod){

  TCanvas * c1 = new TCanvas("c", "c", 800, 800);
  c1->Divide(2,2);

  TString baseName = "cs137";

  TH1F* e0;
  TH1F* e1;
  TH1F* e2;
  TH1F* e3;

  TString fileName = baseName + ".root";
  TFile *rF = new TFile(fileName);

  TString e0Name = baseName + Form("_%i", 4*mod);
  rF->GetObject(e0Name, e0);
  c1->cd(3); // lower left Clover leaf / GRETINA crystal 1 (type B)
  e0->GetXaxis()->SetRangeUser(0,800);
  e0->Draw();

  TString e1Name = baseName + Form("_%i", 4*mod+1);
  rF->GetObject(e1Name, e1);
  c1->cd(4); // lower right Clover leaf / GRETINA crystal 2 (type A)
  e1->GetXaxis()->SetRangeUser(0,800);
  e1->Draw();

  TString e2Name = baseName + Form("_%i", 4*mod+2);
  rF->GetObject(e2Name, e2);
  c1->cd(2); // upper right Clover leaf / GRETINA crystal 3 (type B)
  e2->GetXaxis()->SetRangeUser(0,800);
  e2->Draw();

  TString e3Name = baseName + Form("_%i", 4*mod+3);
  rF->GetObject(e3Name, e3);
  c1->cd(1); // upper left Clover leaf / GRETINA crystal 4 (type A)
  e3->GetXaxis()->SetRangeUser(0,800);
  e3->Draw();

}
