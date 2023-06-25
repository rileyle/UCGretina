#include "Riostream.h"
#include "iostream"
#include "TRandom.h"
#include "TFile.h"

TH1F *clover;
TH1F *crystal[12];
TH2F *crys_xy[12];
TH2F *crys_xz[12];
TH2F *crys_yz[12];
TH1F *crys_x[12];
TH1F *crys_y[12];
TH1F *crys_z[12];

TH1F *crystal_lo[4];
TH2F *crys_xy_lo[4];
TH2F *crys_xz_lo[4];
TH2F *crys_yz_lo[4];
TH1F *crys_x_lo[4];
TH1F *crys_y_lo[4];
TH1F *crys_z_lo[4];

TH1F *crystal_hi[4];
TH2F *crys_xy_hi[4];
TH2F *crys_xz_hi[4];
TH2F *crys_yz_hi[4];
TH1F *crys_x_hi[4];
TH1F *crys_y_hi[4];
TH1F *crys_z_hi[4];

Float_t ecalPar0, ecalPar1;
Float_t sigmaPar0[136] = {0.0};
Float_t sigmaPar1[136] = {0.0};
Float_t sigmaPar2[136] = {0.0};
Float_t sigmaPar3[136] = {0.0};

Float_t threshPar0[136] = {0.0};
Float_t threshPar1[136] = {0.0};

// Read beta, target offset, and resolution parameters
void loadGretina(TString fileName) {

  char line[1000];
  Float_t fBuffer[4];
  TString sBuffer[8];

  ifstream fp;
  fp.open(fileName);

  // Energy calibration parameters
  fp >> ecalPar0  >> ecalPar1;
  fp.getline(line,1000);  // Advance to next line.

  cout << "Energy calibration parameters: "
       << ecalPar0 << ", " << ecalPar1 << endl;

  // Threshold and resolution parameters
  Int_t i;
  while( fp >> i >> threshPar0[i] >> threshPar1[i]
	 >> sigmaPar0[i] >> sigmaPar1[i] 
	 >> sigmaPar2[i] >> sigmaPar3[i] ){
    cout << "Crystal " << i << " threshold parameters: "
	 << threshPar0[i] << ", " << threshPar1[i] << endl;
    cout << "           resolution parameters: "
	 << sigmaPar0[i] << ", " << sigmaPar1[i] 
	 << sigmaPar2[i] << ", " << sigmaPar3[i] 
	 << endl;
    fp.getline(line,1000);  // Advance to next line.
  }
  if(i == 0)
    cout << "No thresholds set!" << endl;

  fp.close();

}

// Simulate intrinsic detector energy resolution and energy (de)calibration
Float_t measuredE(Float_t E, Int_t det){

  E = ecalPar0 + E*ecalPar1;

  Float_t sigma = sigmaPar0[det]
    + sigmaPar1[det]*E
    + sigmaPar2[det]*E*E
    + sigmaPar3[det]*E*E*E;

  return E + gRandom->Gaus(0,sigma);
}

// Read simulation data and sort into histograms
void loadSim(TString fileName) {

  char line[1000];
  Float_t fBuffer[4];
  Int_t   iBuffer[2];
  TString sBuffer;

  Float_t Edep[1000]  = {0.0};
  Float_t mEdep[1000] = {0.0};
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

  TH1F* sim = new TH1F("clover","clover", nChannels, float(lo), float(hi));

  Int_t nHoles = 3;
  Int_t hole[] = { 31, 32, 33 };
  Int_t holeNum[136];

  Int_t nCrystals = nHoles*4;

  Int_t crystalID[136];
  Int_t crystalNum[136];
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
    TString crystalName = "crys";
    TString crystalLabel;
    crystalLabel.Form("_%d", crystalID[i]);
    crystalName += crystalLabel;

    crystal[i] = new TH1F(crystalName,"", nChannels, float(lo), float(hi));

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
    crys_x[i] = new TH1F(xName,"", 1000, -50., 50.);
    TString yName = crystalName.Copy();
    yName += "_y";
    crys_y[i] = new TH1F(yName,"", 1000, -50., 50.);
    TString zName = crystalName.Copy();
    zName += "_z";
    crys_z[i] = new TH1F(zName,"", 1000, 0., 100.);

  }

  // Clover coincidence spectra
  for(Int_t i = 0; i < 4; i++){

    TString crystalName = "crys";
    TString crystalLabel;
    crystalLabel.Form("_%d", crystalID[i]);
    crystalName += crystalLabel;

    TString cName = crystalName.Copy();
    cName += "_lo";
    crystal_lo[i] = new TH1F(cName,"", nChannels, float(lo), float(hi));
    cName = crystalName.Copy();
    cName += "_hi";
    crystal_hi[i] = new TH1F(cName,"", nChannels, float(lo), float(hi));

    TString xyName = crystalName.Copy();
    xyName += "_xy_lo";
    crys_xy_lo[i] = new TH2F(xyName,"", 100, -50., 50., 100, -50., 50.);
    xyName = crystalName.Copy();
    xyName += "_xy_hi";
    crys_xy_hi[i] = new TH2F(xyName,"", 100, -50., 50., 100, -50., 50.);

    TString xzName = crystalName.Copy();
    xzName += "_xz_lo";
    crys_xz_lo[i] = new TH2F(xzName,"", 100, -50., 50., 100, 0., 100.);
    xzName = crystalName.Copy();
    xzName += "_xz_hi";
    crys_xz_hi[i] = new TH2F(xzName,"", 100, -50., 50., 100, 0., 100.);

    TString yzName = crystalName.Copy();
    yzName += "_yz_lo";
    crys_yz_lo[i] = new TH2F(yzName,"", 100, -50., 50., 100, 0., 100.);
    yzName = crystalName.Copy();
    yzName += "_yz_hi";
    crys_yz_hi[i] = new TH2F(yzName,"", 100, -50., 50., 100, 0., 100.);

    TString xName = crystalName.Copy();
    xName += "_x_lo";
    crys_x_lo[i] = new TH1F(xName,"", 1000, -50., 50.);
    xName = crystalName.Copy();
    xName += "_x_hi";
    crys_x_hi[i] = new TH1F(xName,"", 1000, -50., 50.);

    TString yName = crystalName.Copy();
    yName += "_y_lo";
    crys_y_lo[i] = new TH1F(yName,"", 1000, -50., 50.);
    yName = crystalName.Copy();
    yName += "_y_hi";
    crys_y_hi[i] = new TH1F(yName,"", 1000, -50., 50.);

    TString zName = crystalName.Copy();
    zName += "_z_lo";
    crys_z_lo[i] = new TH1F(zName,"", 1000, 0., 100.);
    zName = crystalName.Copy();
    zName += "_z_hi";
    crys_z_hi[i] = new TH1F(zName,"", 1000, 0., 100.);

  }

  for(Int_t i = 0; i < nHoles; i++){

    holeNum[4*hole[i]] = i;
    holeNum[4*hole[i]+1] = i;
    holeNum[4*hole[i]+2] = i;
    holeNum[4*hole[i]+3] = i;

  }
 
  Int_t nGamma = 0;
  Int_t nBuffer = 0;
  Int_t nReaction = 0;
  Int_t nEvent;
  Int_t nHits;
  Int_t fullEnergy;
  Float_t Egamma;

  while( fp >> sBuffer ){

    nBuffer++;
    if(nBuffer%100000 == 0) {
      cout << "Processed " << nBuffer << " events ...\r";
      flush(cout);
    }

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
	if(!(Edep[detNum[k]] > 0.)){
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
    Bool_t cloverLoHit = false;
    Bool_t cloverHiHit = false;
    if(ElabAB > 0) {

      // Simulate intrinsic resolution and low-energy thresholds.
      // Check for Clover hits and target GRETINA crystal hits.
      for(Int_t i = 0; i < nHits; i++) {
	//Simulate thresholds
	if( gRandom->Rndm()
	    < (1.0 + tanh((Edep[detNum[i]]-threshPar0[detNum[i]])/threshPar1[detNum[i]]))/2.0 )
	  mEdep[detNum[i]] = measuredE(Edep[detNum[i]], detNum[i]);
	else
	  mEdep[detNum[i]] = 0.;
	
	if( mEdep[detNum[i]] > 0. ) {
	  if( detNum[i] > 127 ) {
	    if( detNum[i]%4 < 2)         // Bottom-row clover crystal
	      cloverLoHit = true;
	    else                         // Top-row clover crystal
	      cloverHiHit = true;
	  }
	}
      }

      Float_t cloverTotalE = 0.;
      for(Int_t i = 0; i < nHits; i++){
	Float_t E =  Edep[detNum[i]];
	Float_t mE = mEdep[detNum[i]];
	if( mE > 0.0 ){

	  // Cover addback
	  if(detNum[i]>127 && mE > 0.)
	    cloverTotalE += mE;

	  // No addback
	  crystal[crystalNum[detNum[i]]]->Fill( mE );
	  //	    cout << "(" << xHit[detNum[i]] << ", " 
	  //		 << yHit[detNum[i]] << ")" << endl;


	  // cout << "hit " << i << ", detNum[i] = " << detNum[i]
	  // 	 << ", crystalNum[detNum[i]] = " 
	  // 	 << crystalNum[detNum[i]] << endl;

	  crys_xy[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]],
					       yHit[detNum[i]]);
	  crys_xz[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]],
					       zHit[detNum[i]]);
	  crys_yz[crystalNum[detNum[i]]]->Fill(yHit[detNum[i]],
					       zHit[detNum[i]]);

	  crys_x[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]]);
	  crys_y[crystalNum[detNum[i]]]->Fill(yHit[detNum[i]]);
	  crys_z[crystalNum[detNum[i]]]->Fill(zHit[detNum[i]]);

	  // Register clover conicidences
	  if(cloverLoHit && detNum[i] < 128){
	    crystal_lo[crystalNum[detNum[i]]]->Fill( mE );
	    crys_xy_lo[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]],
						    yHit[detNum[i]]);
	    crys_xz_lo[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]],
						    zHit[detNum[i]]);
	    crys_yz_lo[crystalNum[detNum[i]]]->Fill(yHit[detNum[i]],
						    zHit[detNum[i]]);

	    crys_x_lo[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]]);
	    crys_y_lo[crystalNum[detNum[i]]]->Fill(yHit[detNum[i]]);
	    crys_z_lo[crystalNum[detNum[i]]]->Fill(zHit[detNum[i]]);
	  }

	  if(cloverHiHit && detNum[i] < 128){
	    crystal_hi[crystalNum[detNum[i]]]->Fill( mE );
	    crys_xy_hi[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]],
						    yHit[detNum[i]]);
	    crys_xz_hi[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]],
						    zHit[detNum[i]]);
	    crys_yz_hi[crystalNum[detNum[i]]]->Fill(yHit[detNum[i]],
						    zHit[detNum[i]]);

	    crys_x_hi[crystalNum[detNum[i]]]->Fill(xHit[detNum[i]]);
	    crys_y_hi[crystalNum[detNum[i]]]->Fill(yHit[detNum[i]]);
	    crys_z_hi[crystalNum[detNum[i]]]->Fill(zHit[detNum[i]]);
	  }

	}

	nDets++;
	Edep[detNum[i]] = 0.0; // Very important to zero after use.
	xHit[detNum[i]] = 0.;
	yHit[detNum[i]] = 0.;
	zHit[detNum[i]] = 0.;

      }
      if(cloverTotalE > 0.)
	sim->Fill( cloverTotalE );
      nGamma++;
    }
  }
  
  fp.close();

  cout << "\r... sorted " 
       << nBuffer    << " buffers and "
       << nGamma     << " gamma events.\n" << endl;

  return;

};

void ScanSort(TString simFile){

  cout << "\n... loading resolution parameters ...\n" << endl;
  loadGretina("ScanSort.inp");

  TString rootFile = simFile.Copy();
  rootFile.Replace(rootFile.Index(".out",4),4,".root",5);
  TFile *rF = new TFile(rootFile,"recreate");

  TString sourceType = simFile.Copy();
  sourceType.Replace(sourceType.Index(".out",4),4,"",0);

  cout << "... sorting events in " << simFile << " ...\n" << endl;

  loadSim(simFile);

  cout << "... writing spectra to " << rootFile << " ...\n" << endl;

  rF->Write();

  return;

}
