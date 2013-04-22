#include "Riostream.h"
#include "iostream"

// Set /IonPrint/Track_Set in the macro file, and send stdout to a log file. 
// makeTree("filename") crawls the log file and extracts the ion data.

void makeTree(TString fileName) {

  struct ion_s { 

    Float_t KE;
    Float_t betaf;
    Float_t betar;
    Float_t theta;
    Float_t phi;
    Float_t x;
    Float_t y;
    Float_t z;

  };

  ion_s ion;

  TTree *tree = new TTree("tree","Ion Simulation Data");
  tree->Branch("ion",&ion.KE,
	       "KE/F:betaf/F:betar/F:theta/F:phi/F:x/F:y/F:z/F");

  char line[1000];
  char hitMarker[] = "-------->Hits";
  Int_t   iBuffer[16];
  Float_t fBuffer[16];
  TString sBuffer[16];

  ifstream fp;
  fp.open(fileName);

  Int_t nEvt = 0;
  while(fp >> sBuffer[0])
    {
      // Find the next event.
      if( sBuffer[0] == hitMarker ) {

 	nEvt++;
	if(nEvt%1000 == 0) {
	  cout << "Event " << nEvt << endl;
	}

	fp >> sBuffer[1] >> sBuffer[2] >> sBuffer[3] 
	 >> sBuffer[4] >> sBuffer[5] >> sBuffer[6] >> iBuffer[0];

	int nHits = iBuffer[0];

	fp.getline(line,1000); // Advance to the next line.
	fp.getline(line,1000); // Skip titles.
	//	if(!(fp.good())) break;

	// Process the hit descriptions
	for(int i=0; i < nHits; i++) {
	  fp >> iBuffer[0] >> iBuffer[1] >> sBuffer[0] >> fBuffer[0] 
             >> fBuffer[1] >> fBuffer[2] >> fBuffer[3] >> fBuffer[4] 
             >> fBuffer[5] >> fBuffer[6];

	  //cout << iBuffer[0] << iBuffer[1] << sBuffer[0] << fBuffer[0] 
	  //     << fBuffer[1] << fBuffer[2] << fBuffer[3] << fBuffer[4] 
	  //     << fBuffer[5] << fBuffer[6] << endl;
	  
	  if(iBuffer[1] == 5) {
	    ion.betar = fBuffer[1];
	    ion.x = fBuffer[4];
	    ion.y = fBuffer[5];
	    ion.z = fBuffer[6];
	  }

	  if(i==nHits-1) {
	    ion.KE    = fBuffer[0];
	    ion.betaf = fBuffer[1];
	    ion.theta = fBuffer[2];
	    ion.phi   = fBuffer[3];
	  }

	  fp.getline(line,1000);  // Move on to the next line.

	}

	//cout << ion.KE << " "  << ion.beta << " " << ion.theta << " "
	//     << ion.phi << " " << ion.x << " " << ion.y << " " << ion.z
        //     << endl;

	tree->Fill();

	fp.getline(line,1000);    // Discard the blank line between events.

      }

    }
  fp.close();

  tree->Print();
  //f1->Write();
  
  tree->StartViewer();
  return;

};
