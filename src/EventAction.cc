#include "EventAction.hh"
#include "RunAction.hh"

EventAction::EventAction()
{ 
  ionCollectionID=-1;
  gammaCollectionID=-1;
  positionRes = 0.*mm;
  crmatFileName = "";
  gretinaCoords = false;
}


EventAction::~EventAction()
{
  ;
}

void EventAction::BeginOfEventAction(const G4Event*)
{
 
  G4SDManager * SDman = G4SDManager::GetSDMpointer();

  if(gammaCollectionID<0||ionCollectionID<0)
    {
      gammaCollectionID=SDman->GetCollectionID("gammaCollection");
      ionCollectionID=SDman->GetCollectionID("ionCollection");
    }

  // G4cout<<"+++++ Begin of event "<<evt->GetEventID()<<G4endl;

}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  // G4cout<<"+++++ End of event "<<evt->GetEventID()<<G4endl;

  G4int event_id=evt->GetEventID();

  if(event_id%10000==0) {
    G4cout<<" Number of processed events "<<event_id<<"\r"<<std::flush;
  }
 
  // Write event information to the output file.
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  if(HCE) {

    TrackerGammaHitsCollection* gammaCollection = (TrackerGammaHitsCollection*)(HCE->GetHC(gammaCollectionID));

    G4int Nhits = gammaCollection->entries();

    if(Nhits>0) {

      // Consolidate total each gamma-ray interaction with its children 
      // based on proximity. 

      G4int detNum[1000];
      G4double measuredEdep[1000];
      G4double measuredX[1000];
      G4double measuredY[1000];
      G4double measuredZ[1000];
      G4int NCons[1000];
      G4double positionRes2 = positionRes*positionRes;

      // Initialize the first "measured" interaction point.

      detNum[0] = (*gammaCollection)[0]->GetDetNumb();
      measuredEdep[0] = (*gammaCollection)[0]->GetEdep()/keV;
      measuredX[0] = (*gammaCollection)[0]->GetPos().getX()/mm;
      measuredY[0] = (*gammaCollection)[0]->GetPos().getY()/mm;
      measuredZ[0] = (*gammaCollection)[0]->GetPos().getZ()/mm;
      NCons[0] = 1;
      G4int NMeasured = 1;

      G4double totalEdep = (*gammaCollection)[0]->GetEdep()/keV;

      for(G4int i = 1; i < Nhits; i++){

	totalEdep += (*gammaCollection)[i]->GetEdep()/keV;

	// Compare hit i with existing "measured" interaction points.
	// Only consolidate within single crystals.

	G4double x = (*gammaCollection)[i]->GetPos().getX()/mm;
	G4double y = (*gammaCollection)[i]->GetPos().getY()/mm;
	G4double z = (*gammaCollection)[i]->GetPos().getZ()/mm;
	G4double e = (*gammaCollection)[i]->GetEdep()/keV;

	G4bool consolidated = false;

	for(G4int j = 0; j < NMeasured; j++){

	  G4double prox2 = (x - measuredX[j])*(x - measuredX[j]) + (y - measuredY[j])*(y - measuredY[j]) + (z - measuredZ[j])*(z - measuredZ[j]);

	  if( prox2 < positionRes2                                   // proximal
	     && (*gammaCollection)[i]->GetDetNumb() == detNum[j] // same crystal
	     && !consolidated ){                          // not already counted

	    // Energy-weighted average
	    measuredX[j] = (measuredEdep[j]*measuredX[j] + e*x)/(measuredEdep[j] + e);
	    measuredY[j] = (measuredEdep[j]*measuredY[j] + e*y)/(measuredEdep[j] + e);
	    measuredZ[j] = (measuredEdep[j]*measuredZ[j] + e*z)/(measuredEdep[j] + e);
	    measuredEdep[j] += e;

	    NCons[j]++;
	    consolidated = true;

	  }

	}

	// If hit i cannot be consolidated with an existing gamma-ray 
	// interaction point, add a new one. 

	if(!consolidated){

	  detNum[NMeasured] = (*gammaCollection)[i]->GetDetNumb();
	  measuredEdep[NMeasured] = e;
	  measuredX[NMeasured] = x;
	  measuredY[NMeasured] = y;
	  measuredZ[NMeasured] = z;

	  NCons[NMeasured] = 1;
	  NMeasured++;

	}

      }

      // Consolidate the "measured" gamma-ray interaction points.
      // (In some cases, the tracks of secondary particles blur the
      // distinction between initial gamma-ray hits.) 
      G4int NGammaHits = NMeasured;
      for(G4int i = 0; i < NMeasured; i++){
	for(G4int j = i+1; j < NMeasured; j++){

	  if( ( (measuredX[i] - measuredX[j])*(measuredX[i] - measuredX[j])
		   + (measuredY[i] - measuredY[j])*(measuredY[i] - measuredY[j])
		   + (measuredZ[i] - measuredZ[j])*(measuredZ[i] - measuredZ[j])
		   < positionRes2 )                                  // proximal
	      && detNum[i] == detNum[j]                          // same crystal
	      &&  (NCons[i] > 0 && NCons[j] > 0) ){  // not already consolidated

	    // Energy-weighted average
	    measuredX[i] = (measuredEdep[i]*measuredX[i] + measuredEdep[j]*measuredX[j])/(measuredEdep[i]+measuredEdep[j]);
	    measuredY[i] = (measuredEdep[i]*measuredY[i] + measuredEdep[j]*measuredY[j])/(measuredEdep[i]+measuredEdep[j]);
	    measuredZ[i] = (measuredEdep[i]*measuredZ[i] + measuredEdep[j]*measuredZ[j])/(measuredEdep[i]+measuredEdep[j]);
	    measuredEdep[i] += measuredEdep[j];

	    NCons[j] = -1;
	    NGammaHits--;

	  }
	  
	}

      }

      // Write the event with consolidated interaction points to the
      // output file. 
      evfile << "$   " << NGammaHits << "   " 
	     << totalEdep << "   "
	     << event_id
	     << endl;

      for(int i=0; i<NMeasured; i++) {

	if(NCons[i] > 0){

	  G4double x, y, z;
	  if(gretinaCoords){
	    x = -measuredY[i];      // If the user specified GRETINA 
	    y = measuredX[i];       // coordinates, rotate pi/2 about z.
	    measuredX[i] = x;
	    measuredY[i] = y;
	  } else {
	    x = measuredX[i];
	    y = measuredY[i];
	  }
	  z = measuredZ[i];

	  // If the user supplied a crmat (crystal -> world coords), 
	  // invert the transformation (world -> crystal).
	  if(crmatFileName != "") {
	    G4int h = (G4int)detNum[i]/4;  // (Hole-1)
	    G4int c = detNum[i]%4;         // Crystal

	    // Reverse transformation: first subtract the translation ...
	    x -= crmat[h][c][0][3];
	    y -= crmat[h][c][1][3];
	    z -= crmat[h][c][2][3];
	    // ... then apply the inverse (transpose) rotation matrix.
	    measuredX[i] = x*crmat[h][c][0][0] 
	                 + y*crmat[h][c][1][0]
                         + z*crmat[h][c][2][0];
	    measuredY[i] = x*crmat[h][c][0][1] 
	                 + y*crmat[h][c][1][1]
                         + z*crmat[h][c][2][1];
	    measuredZ[i] = x*crmat[h][c][0][2] 
	                 + y*crmat[h][c][1][2]
                         + z*crmat[h][c][2][2];
	  }

	  evfile
	    << std::setw(3) << detNum[i]+4  // +4 to match measured data stream
	    << std::fixed << std::setprecision(2) << std::setw(12) << std::right
	    << measuredEdep[i]
	    << std::setw(10) << std::right
	    << measuredX[i]
	    << std::setw(10) << std::right
	    << measuredY[i]
	    << std::setw(10) << std::right
	    << measuredZ[i]
	    << G4endl;

	}

      }

    }

  }

}

// --------------------------------------------------TB
void EventAction::openEvfile()
{
  //LR: Already open if name is set in macro file, so check:
  if (!evfile.is_open()) evfile.open(outFileName.c_str());

  if (evfile == NULL) G4cout<< "evfile ERROR"<<G4endl;
  
  return;
}
//--------------------------------------------------- TB
void EventAction::closeEvfile()
{
	evfile.close();
	return;
}
//----------------------------------------------------TB
void EventAction::SetOutFile(G4String name)
{
	outFileName = name;
	G4cout << "\nOutput file: " << outFileName << G4endl;
	closeEvfile();
	openEvfile();
	return;
}
// --------------------------------------------------
void EventAction::openCrmatFile()
{
  if (!crmatFile.is_open()) crmatFile.open(crmatFileName.c_str());

  if (crmatFile == NULL) G4cout << "crmatFile ERROR" << G4endl;

  return;
}
//---------------------------------------------------
void EventAction::closeCrmatFile()
{
  crmatFile.close();
  return;
}
//---------------------------------------------------
void EventAction::SetCrmatFile(G4String name) { 

  crmatFileName = name;

  openCrmatFile();

  G4cout << "\nUsing crmat from file " << crmatFileName << " to transform interaction points from world to crytsal coordinates." << G4endl;

  G4int hole, crystal;
  for(G4int h = 0; h < 30; h++){
    for(G4int c = 0; c < 4; c++){
      crmatFile >> hole >> crystal;
      //      G4cout << "Hole: " << hole << "  Crystal: " << crystal << G4endl;
      if(hole != h || crystal != c){
	G4cout << "Error reading crmat file " << crmatFileName << G4endl;
	exit(-1);
      }
      for(G4int i = 0; i < 4; i++){
	for(G4int j = 0; j < 4; j++){
	  crmatFile >> crmat[h][c][i][j];
	  if (j == 3)
	    crmat[h][c][i][3] *= cm;  	          // Convert translations to mm.
	  //	  G4cout << crmat[h][c][i][j] << "   ";
	}
	//	G4cout << G4endl;
      }
    }
  }

  closeCrmatFile();

  return;
}
//---------------------------------------------------
void EventAction::SetGretinaCoords(){

  gretinaCoords = true; 

  G4cout << "Writing interaction points in the Gretina coordinate system (x = down, z = beam)." << G4endl;

  return;
}
//---------------------------------------------------
