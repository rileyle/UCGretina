#include "Riostream.h"
#include "iostream"
#include "TRandom.h"
#include "TFile.h"

TH1F *sim;
TH1F *sim_addback;
TH1F *crystal[4];
TH2F *crys_xy[4];
TH2F *crys_xz[4];
TH2F *crys_yz[4];
TH1F *crys_x[4];
TH1F *crys_y[4];
TH1F *crys_z[4];

TH1F *mod_x[32];
TH1F *mod_y[32];
TH1F *mod_z[32];

Float_t sigmaPar0, sigmaPar1, sigmaPar2, sigmaPar3;
Float_t ecalPar0, ecalPar1;
Float_t threshPar0[128] = {0.0};
Float_t threshPar1[128] = {0.0};

// Read beta, target offset, and resolution parameters
void loadGretina(TString fileName) {

  char line[1000];
  Float_t fBuffer[4];
  TString sBuffer[8];

  ifstream fp;
  fp.open(fileName);

  // Resolution parameters
  fp >> sigmaPar0 >> sigmaPar1 >> sigmaPar2 >> sigmaPar3;
  fp.getline(line,1000);  // Advance to next line.

  // Energy calibration parameters
  fp >> ecalPar0  >> ecalPar1;
  fp.getline(line,1000);  // Advance to next line.

  // Threshold parameters
  Int_t i;
  while(fp >> i >> threshPar0[i] >> threshPar1[i]){
    cout << "Crystal " << i << " threshold parameters: "
	 << threshPar0[i] << ", " << threshPar1[i] << endl;
    fp.getline(line,1000);  // Advance to next line.
  }
  if(i == 0)
    cout << "No thresholds set!" << endl;

  fp.close();

}

// Simulate intrinsic detector energy resolution and energy (de)calibration
Float_t measuredE(Float_t E){

  E = E*ecalPar1 + ecalPar0;

  Float_t sigma = sigmaPar0 
    + sigmaPar1*E
    + sigmaPar2*E*E
    + sigmaPar3*E*E*E;

  return E + gRandom->Gaus(0,sigma);
}

// Read simulation data and sort into histograms
void loadSim(TString fileName) {

  char line[1000];
  Float_t fBuffer[4];
  Int_t   iBuffer[2];
  TString sBuffer;

  Float_t Edep[1000] = {0.0};
  Float_t xHit[1000] = {0.0};
  Float_t yHit[1000] = {0.0};
  Float_t zHit[1000] = {0.0};
  Int_t   detNum[200] = {0};

  ifstream fp;
  fp.open(fileName);

  Int_t nChannels     = 4096;
  Int_t keVperChannel = 1;    // Match binning of the spectrum to fit
  Int_t lo = 0;
  Int_t hi = lo+nChannels*keVperChannel;

  TString Name = fileName.Copy();
  Name.Replace(Name.Index(".out",4,0,0),4,"",0);
  sim = new TH1F(Name,"", nChannels, float(lo), float(hi));

  TString NameAB = fileName.Copy();
  NameAB.Replace(NameAB.Index(".out",4,0,0),4,"_addback",8);
  sim_addback = new TH1F(NameAB,"", nChannels, float(lo), float(hi));


  Int_t nHoles = 1;
  Int_t hole[] = { 31 };
  Int_t holeNum[128];

  Int_t nCrystals = nHoles*4;

  Int_t crystalID[128];
  Int_t crystalNum[128];
  for(Int_t i = 0; i < nHoles; i++){
    crystalID[4*i] = 4*hole[i];
    crystalID[4*i+1] = 4*hole[i]+1;
    crystalID[4*i+2] = 4*hole[i]+2;
    crystalID[4*i+3] = 4*hole[i]+3;
  }

  for(Int_t i = 0; i < nCrystals; i++){
    // Load detector ID lookup table
    crystalNum[crystalID[i]] = i;

    // Initialize crystal spectra
    TString crystalName = Name.Copy();
    TString crystalLabel;
    //    crystalLabel.Form("_%02d", i);
    crystalLabel.Form("_%d", crystalID[i]);
    crystalName += crystalLabel;
    crystal[i] = new TH1F(crystalName,"", nChannels, float(lo), float(hi));
    //    crystal[i]->Sumw2();
    TString xyName = crystalName.Copy();
    xyName += "_xy";
    crys_xy[i] = new TH2F(xyName,"", 100, -50., 50., 100, -50., 50.);

    TString xzName = crystalName.Copy();
    xzName += "_xz";
    crys_xz[i] = new TH2F(xzName,"", 100, -50., 50., 100, 0., 100.);

    TString yzName = crystalName.Copy();
    yzName += "_yz";
    crys_yz[i] = new TH2F(yzName,"", 100, -50., 50., 100, 0., 100.);

    TString xName = crystalName.Copy();
    xName += "_x";
    crys_x[i] = new TH1F(xName,"", 1000, -500., 500.);
    TString yName = crystalName.Copy();
    yName += "_y";
    crys_y[i] = new TH1F(yName,"", 1000, -500., 500.);
    TString zName = crystalName.Copy();
    zName += "_z";
    crys_z[i] = new TH1F(zName,"", 1000, -500., 500.);

  }

  for(Int_t i = 0; i < nHoles; i++){

    holeNum[4*hole[i]] = i;
    holeNum[4*hole[i]+1] = i;
    holeNum[4*hole[i]+2] = i;
    holeNum[4*hole[i]+3] = i;

    TString modName = Name.Copy();
    TString modLabel;
    modLabel.Form("_mod%02d", hole[i]);
    modName += modLabel;

    TString mxName = modName.Copy();
    mxName += "_x";
    mod_x[i] = new TH1F(mxName,"", 1000, -500., 500.);

    TString myName = modName.Copy();
    myName += "_y";
    mod_y[i] = new TH1F(myName,"", 1000, -500., 500.);

    TString mzName = modName.Copy();
    mzName += "_z";
    mod_z[i] = new TH1F(mzName,"", 1000, -500., 500.);

  }
 
  Int_t nGamma = 0;
  Int_t nReaction = 0;
  Int_t nPhotopeak = 0;
  Int_t nGamma, nEvent;
  Int_t nHits;
  Int_t fullEnergy;
  Float_t Egamma;

  while( fp >> sBuffer ){

    Int_t Ndecomp, eventID;
    if( sBuffer == "D" )
      fp >> Ndecomp >> eventID;
    else
      continue;

    //    cout << "Ndecomp = " << Ndecomp << "  eventID = " << eventID << endl;

    Float_t ElabAB = 0;
    Float_t x,   y, z;
    Float_t xs, ys, zs;
    Float_t EdepMax = 0;

    Float_t Edeposited;

    // Load the Edep array: energy deposited, indexed by detector number
    //    and detNum array: detector number, indexed by hit number
    Int_t k = 0;
    for(Int_t i = 0; i < Ndecomp; i++){

      Int_t Nip;

      fp >> sBuffer >> iBuffer[0] >> iBuffer[1];

      detNum[k] = iBuffer[0];
      Nip = iBuffer[1];

      //      cout << sBuffer 
      //	   << " detNum[" << k << "] = " << detNum[k] 
      //	   << "  Nip = " << Nip << endl;

      for(Int_t j = 0; j < Nip; j++){

	fp >> iBuffer[0] >> Edeposited >> x >> y >> z;

	//	cout << "Edeposited = " << Edeposited 
	//	     << "  x = " << x
	//	     << "  y = " << y
	//	     << "  z = " << z
	//	     << endl;

	// Addback
	ElabAB += Edeposited;
	//    Set the hit position = position of largest energy deposit.
	if(Edeposited > EdepMax){
	  EdepMax = Edeposited;
	  xs = x;
	  ys = y;
	  zs = z;
	}

	// No addback
	//   Hit position = first hit in crystal
	if(!Edep[detNum[k]] > 0.){
	  xHit[detNum[k]] = x;
	  yHit[detNum[k]] = y;
	  zHit[detNum[k]] = z;
	}
	Edep[detNum[k]] += Edeposited;

      }

      k++;

    }

    nHits = k;

    Int_t nDets = 0;

    if(ElabAB > 0) {
      ElabAB = 0; // (We'll recalculate it with simulated thresholds.)
      for(Int_t i = 0; i < nHits; i++){
	Float_t E = Edep[detNum[i]];
	if( E > 0.0 ){

	  // Simulate thresholds
	  if( gRandom->Rndm()
	      < (1.0 + tanh((E-threshPar0[detNum[i]])/threshPar1[detNum[i]]))/2.0 ){

	    // No addback
	    sim->Fill( measuredE(E) );
	    // No addback
	    crystal[crystalNum[detNum[i]]]->Fill( measuredE(E) );
	    //	    cout << "(" << xHit[detNum[i]] << ", " 
	    //		 << yHit[detNum[i]] << ")" << endl;

	    crys_xy[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]],
						 yHit[detNum[i]]);

	    crys_xz[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]],
						 zHit[detNum[i]]);

	    crys_yz[crystalNum[detNum[i]]]->Fill(yHit[detNum[i]],
						 zHit[detNum[i]]);

	    Float_t xH = xHit[detNum[i]];
	    Float_t yH = yHit[detNum[i]];
	    Float_t zH = zHit[detNum[i]];

	    crys_x[crystalNum[detNum[i]]]->Fill(xH);
	    crys_y[crystalNum[detNum[i]]]->Fill(yH);
	    crys_z[crystalNum[detNum[i]]]->Fill(zH);

	    mod_x[holeNum[detNum[i]]]->Fill(xH);
	    mod_y[holeNum[detNum[i]]]->Fill(yH);
	    mod_z[holeNum[detNum[i]]]->Fill(zH);

	    // Addback			
	    ElabAB += E;

	  }

	  nDets++;
	  Edep[detNum[i]] = 0.0; // Very important to zero after use.
	  xHit[detNum[i]] = 0.;
	  yHit[detNum[i]] = 0.;
	  zHit[detNum[i]] = 0.;

	}
      }

      // Addback
      sim_addback->Fill( measuredE(ElabAB) );

      // Photopeak events
      if( abs(Egamma - ElabAB) < 0.02 )
	nPhotopeak++;
    }

    nGamma++;
    if(nGamma%100000 == 0) {
      cout << "\r... sorted " << nGamma << " gamma rays ...";
      flush(cout);
    }

  }

  fp.close();

  cout << "\r... sorted " 
       << nGamma     << " gamma events, and " 
       << nPhotopeak << " photopeak events ...\n" << endl;

  return;

};

void ScanSort(TString simFile){

  cout << "\n... loading resolution parameters ...\n" << endl;
  loadGretina("ScanSort.inp");

  TString rootFile = simFile.Copy();
  rootFile.Replace(rootFile.Index(".out",4,0,0),4,".root",5);
  TFile *rF = new TFile(rootFile,"recreate");

  TString sourceType = simFile.Copy();
  sourceType.Replace(sourceType.Index(".out",4,0,0),4,"",0);

  cout << "... sorting events in " << simFile << " ...\n" << endl;

  loadSim(simFile);

  cout << "... writing Spectra to " << rootFile << " ...\n" << endl;

  rF->Write();

  return;

}
