#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Timer.hh"
extern G4Timer Timerintern;

EventAction::EventAction()
{ 
  ionCollectionID=-1;
  gammaCollectionID=-1;
  hitTolerance = 0.00001*mm;
  packingRes = 0.*mm;
  S800KE = 1.0;
  allS800 = false;
  outFileName = "";
  outDetsOnly = false;
  evOut = false;
  mode2FileName = "";
  mode2Out = false;
  crmatFileName = "";
  crystalXforms = false;
  gretinaCoords = false;
  print = false;
  evt = NULL;
  fisInBeam = false;
  timeSort = false;
  Timerintern.Start();
  timerCount = 0;
  eventsPerSecond = 0;
  everyNevents = 1000;
  posRes = 0.;
  threshE = 0.;
  threshDE = 0.001*keV;
}


EventAction::~EventAction()
{
  ;
}

void EventAction::BeginOfEventAction(const G4Event* ev)
{
  evt = ev;

  // Add a G4UserEventInformation object to store event info
  EventInformation* eventInfo = new EventInformation;
  G4EventManager::
    GetEventManager()->SetUserInformation(eventInfo);

  G4SDManager * SDman = G4SDManager::GetSDMpointer();

  if(gammaCollectionID<0||ionCollectionID<0)
    {
      gammaCollectionID=SDman->GetCollectionID("gammaCollection");
      ionCollectionID=SDman->GetCollectionID("ionCollection");
    }

  // For event filter
  eventInfo->SetWriteEvent(false);
  
  // G4cout<<"+++++ Begin of event "<<evt->GetEventID()<<G4endl;

}


void EventAction::EndOfEventAction(const G4Event* ev)
{
  evt = ev;

  G4int event_id=evt->GetEventID();

  if(event_id%everyNevents == 0 && event_id > 0) {

    std::ios::fmtflags f( G4cout.flags() );
    G4int prec = G4cout.precision();

    Timerintern.Stop();
    timerCount++;
    eventsPerSecond += 
      ((double)everyNevents/Timerintern.GetRealElapsed() 
       - eventsPerSecond)/timerCount;
    G4cout << std::fixed << std::setprecision(0) << std::setw(3) 
	   << std::setfill(' ')
	   << (float)event_id/NTotalEvents*100 << " %   "
	   << eventsPerSecond << " events/s ";

    G4double hours, minutes, seconds;
    G4double time = (float)(NTotalEvents - event_id)/eventsPerSecond;
    hours = floor(time/3600.0);
    if(hours>0){
      G4cout << std::setprecision(0) << std::setw(2) 
	     << hours << ":";
      G4cout << std::setfill('0');
    } else {
      G4cout << std::setfill(' ');
    }
    minutes = floor(fmod(time,3600.0)/60.0);
    if(minutes>0){
      G4cout << std::setprecision(0) << std::setw(2) << minutes << ":";
      G4cout << std::setfill('0');
    } else {
      G4cout << std::setfill(' ');
    }
    seconds = fmod(time,60.0);
    if(seconds>0)
      G4cout << std::setprecision(0) << std::setw(2) << seconds;
    G4cout << std::setfill(' ');
    G4cout << " remaining       "
	   << "\r"<<std::flush;
    Timerintern.Start();

    G4cout.setf( f );
    G4cout.precision( prec );

  }
  
  EventInformation* eventInfo = (EventInformation*)evt->GetUserInformation();

  // Event filter
  if(!eventInfo->WriteEvent())
    return;
  
  // All Mode 2 output from this event gets this timestamp.
  long long int timestamp = (long long int)10000*event_id;

  if(print){
    
    std::ios::fmtflags f( G4cout.flags() );
    G4int prec = G4cout.precision();

    G4cout << "-------->Mode2 data, event " << event_id << G4endl;
    if(fisInBeam)
      G4cout << std::fixed << std::setprecision(4) 
	     << std::setw(12) << std::right
	     << " ata = " << eventInfo->GetATA()
	     << " bta = " << eventInfo->GetBTA()
	     << " dta = " << eventInfo->GetDTA()
	     << " yta = " << eventInfo->GetYTA()
	     << G4endl;
    G4cout << eventInfo->GetNEmittedGammas() << " emitted gamma(s)" << G4endl;
    for(G4int i = 0; i< eventInfo->GetNEmittedGammas(); i++)
      G4cout << std::fixed << std::setprecision(4) 
	     << std::setw(12) << std::right
	     << "energy = " << eventInfo->GetEmittedGammaEnergy(i)
	     << " pos = " << eventInfo->GetEmittedGammaPosX(i)
	     << ", " << eventInfo->GetEmittedGammaPosY(i)
	     << ", " << eventInfo->GetEmittedGammaPosZ(i)
	     << " direction = " << eventInfo->GetEmittedGammaPhi(i)
	     << ", " << eventInfo->GetEmittedGammaTheta(i)
	     << " beta = " << eventInfo->GetBeta(i)
	     << G4endl;

    G4cout.setf( f );
    G4cout.precision( prec );

  }

  // Analyze hits and write event information to the output file.
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  if(HCE) {

    TrackerGammaHitsCollection* gammaCollection 
      = (TrackerGammaHitsCollection*)(HCE->GetHC(gammaCollectionID));

    G4int Nhits = gammaCollection->entries();

    if(Nhits>0) {

      // Packing: consolidate interaction points within segments 
      // based on proximity. 
      G4int trackID[100*MAX_INTPTS];
      G4int detNum[100*MAX_INTPTS];
      G4int segNum[100*MAX_INTPTS];
      G4double measuredEdep[100*MAX_INTPTS];
      G4double segmentEdep[100*MAX_INTPTS];
      G4double measuredX[100*MAX_INTPTS];
      G4double measuredY[100*MAX_INTPTS];
      G4double measuredZ[100*MAX_INTPTS];
      G4double X0[100*MAX_INTPTS];
      G4double Y0[100*MAX_INTPTS];
      G4double Z0[100*MAX_INTPTS];
      G4double globalTime[100*MAX_INTPTS];
      G4int NCons[100*MAX_INTPTS];
      G4double packingRes2 = packingRes*packingRes;

      G4int NMeasured = 0;
      G4double totalEdep = 0;

      // Prevent seg fault in events with very large hit collections.
      // This arises rarely in high-energy muon (and likely other
      // high-energy charged particle) events.
      if( Nhits > 3000 ){
	G4cerr << "Warning: hit collection with " << Nhits
	       << " entries. Processing the first 3000."
	       << " (event " << event_id << ")"
	       << G4endl;
	Nhits = 3000;
      }
      
      for(G4int i = 0; i < Nhits; i++){

	G4double x, y, z;
	if(crystalXforms){
	  G4double xx, yy;
	  xx = (*gammaCollection)[i]->GetPosCrys().getX()/mm;
	  yy = (*gammaCollection)[i]->GetPosCrys().getY()/mm;
	  z = (*gammaCollection)[i]->GetPosCrys().getZ()/mm;
	  G4double ph = -122.93*3.14159/180.;              // Type A
	  if((*gammaCollection)[i]->GetDetNumb()%2 == 0) // Type B
	    ph = -208.25*3.14159/180.;
	  x = xx*cos(ph) - yy*sin(ph);
	  y = xx*sin(ph) + yy*cos(ph);

	} else {
	  x = (*gammaCollection)[i]->GetPos().getX()/mm;
	  y = (*gammaCollection)[i]->GetPos().getY()/mm;
	  z = (*gammaCollection)[i]->GetPos().getZ()/mm;
	}
	G4double en  = (*gammaCollection)[i]->GetEdep()/keV;
	totalEdep += en;

	G4double gt  = (*gammaCollection)[i]->GetGlobalTime()*1.e3;
	
	NCons[i] = -1;
	G4bool processed = false;	

	// Initialize a new interaction point for each gamma-ray hit.
	if((*gammaCollection)[i]->GetParticleID() == "gamma"){

	  // Combine multiple gamma hits at the same position.
	  // (This is rare, but it happens.)
	  if(i > 0 
	     && (x - measuredX[i-1])*(x - measuredX[i-1]) < hitTolerance*hitTolerance
	     && (y - measuredY[i-1])*(y - measuredY[i-1]) < hitTolerance*hitTolerance
	     && (z - measuredZ[i-1])*(z - measuredZ[i-1]) < hitTolerance*hitTolerance){

	    measuredEdep[NMeasured-1] += en; 
	    processed = true;

	  } else {

	    trackID[NMeasured]      = (*gammaCollection)[i]->GetTrackID();
	    detNum[NMeasured]       = (*gammaCollection)[i]->GetDetNumb();
	    segNum[NMeasured]       = (*gammaCollection)[i]->GetSegNumb();

	    // This becomes the total energy deposit associated with this 
	    // interaction.
	    measuredEdep[NMeasured] = en; 

	    // This becomes the barycenter of all energy depositions associated
	    // with this interaction.
	    measuredX[NMeasured]    = x;
	    measuredY[NMeasured]    = y;
	    measuredZ[NMeasured]    = z;
	    trackID[NMeasured]      = (*gammaCollection)[i]->GetTrackID();

	    // Position of the initial interaction. We use position to identify
	    // the tracks produced by this interaction.
	    X0[NMeasured] = (*gammaCollection)[i]->GetPos().getX()/mm;
	    Y0[NMeasured] = (*gammaCollection)[i]->GetPos().getY()/mm;
	    Z0[NMeasured] = (*gammaCollection)[i]->GetPos().getZ()/mm;

	    globalTime[NMeasured] = (*gammaCollection)[i]->GetGlobalTime()*1.e3;
	    
	    NCons[NMeasured] = 1;
	    NMeasured++;
	    processed = true;
	  }

	// Combine secondary-particle hits with their parent interaction points.
	} else {

	  // Compare hit i with existing "measured" interaction points.
	  for(G4int j = 0; j < NMeasured; j++){

	    G4double x0 = (*gammaCollection)[i]->GetTrackOrigin().getX()/mm;
	    G4double y0 = (*gammaCollection)[i]->GetTrackOrigin().getY()/mm;
	    G4double z0 = (*gammaCollection)[i]->GetTrackOrigin().getZ()/mm;

	    // G4cout << "(*gammaCollection)["<< i <<"]->GetParentTrackID() = "
	    // 	   << (*gammaCollection)[i]->GetParentTrackID()
	    // 	   << "   trackID[j] = " << trackID[j]
	    // 	   << "   (x0 - X0[" << j << "]) = " << (x0 - X0[j])
	    // 	   << "   (y0 - Y0[" << j << "]) = " << (y0 - Y0[j])
	    // 	   << "   (z0 - Z0[" << j << "]) = " << (z0 - Z0[j])
	    // 	   << "   (*gammaCollection)[" << i << "]->GetDetNumb()" 
	    // 	   << (*gammaCollection)[i]->GetDetNumb()
	    // 	   << "   detNum[" << j << "] = " << detNum[j] << G4endl;

	    if( (*gammaCollection)[i]->GetParentTrackID() == trackID[j]  // correct parent
		&& (x0 - X0[j])*(x0 - X0[j]) < hitTolerance*hitTolerance
		&& (y0 - Y0[j])*(y0 - Y0[j]) < hitTolerance*hitTolerance
		&& (z0 - Z0[j])*(z0 - Z0[j]) < hitTolerance*hitTolerance // correct interaction point
		&& (*gammaCollection)[i]->GetDetNumb() == detNum[j]){        // same crystal

	      // Energy-weighted average position (barycenter)
	      if(measuredEdep[j] == 0 && en == 0){
		// G4cout << "    Both energies are zero. Taking the average." << G4endl;
		measuredX[j]  = (measuredX[j] + x)/2.;
		measuredY[j]  = (measuredY[j] + y)/2.;
		measuredZ[j]  = (measuredZ[j] + z)/2.;
		
		// Assign the earliest global time of the first raw hit in this IP
		globalTime[j] = std::min(gt, globalTime[j]);
	      } else {
		// G4cout << "    Calculating a weighted average." << G4endl;
		measuredX[j] = (measuredEdep[j]*measuredX[j] + en*x)/(measuredEdep[j] + en);
		measuredY[j] = (measuredEdep[j]*measuredY[j] + en*y)/(measuredEdep[j] + en);
		measuredZ[j] = (measuredEdep[j]*measuredZ[j] + en*z)/(measuredEdep[j] + en);
		measuredEdep[j] += en;

		// Assign the earliest global time of the first raw hit in this IP
		globalTime[j] = std::min(gt, globalTime[j]);
	      }
	      
	      NCons[j]++;
	      processed = true;
	      
	    }
	  }
	}

	// If hit i is not a gamma-ray hit and cannot be consolidated with 
	// an existing gamma-ray interaction point, it's a positron or an
	// electron multiple-scattering event that isn't associated with 
	// energy deposition in the sensitive volume of GRETINA by a gamma 
	// ray. (The gamma-ray interaction happened in dead material.)
	// We'll initialize a new interaction point and treat it as a 
	// gamma-ray interaction.
	if(!processed){

	  trackID[NMeasured]      = (*gammaCollection)[i]->GetTrackID();
	  detNum[NMeasured]       = (*gammaCollection)[i]->GetDetNumb();
	  segNum[NMeasured]       = (*gammaCollection)[i]->GetSegNumb();
	  measuredEdep[NMeasured] = en;
	  measuredX[NMeasured]    = x;
	  measuredY[NMeasured]    = y;
	  measuredZ[NMeasured]    = z;

	  // This is not a gamma ray. We need to trick its siblings into 
	  // treating it as the parent gamma.
	  trackID[NMeasured]      = (*gammaCollection)[i]->GetParentTrackID();
	  X0[NMeasured] = (*gammaCollection)[i]->GetTrackOrigin().getX()/mm;
	  Y0[NMeasured] = (*gammaCollection)[i]->GetTrackOrigin().getY()/mm;
	  Z0[NMeasured] = (*gammaCollection)[i]->GetTrackOrigin().getZ()/mm;

	  globalTime[NMeasured] = (*gammaCollection)[i]->GetGlobalTime()*1.e3;

	  NCons[NMeasured] = 1;
	  NMeasured++;
	  processed = true;

	}

	if(!processed)
	  G4cerr << "Warning: Could not find a home for hit " << i
		 << " of event " << event_id 
		 << G4endl;

	if(NMeasured >= 100*MAX_INTPTS){
	  G4cout << "Error: too many decomposed hits. Increase hit processing array dimension."
		 << G4endl;
	  exit(EXIT_FAILURE);
	}
	
      }

      // Packing: Consolidate the "measured" gamma-ray interaction points
      // within a single segment that are closer than the PackingRes 
      // parameter.
      G4int NGammaHits = NMeasured;
      for(G4int i = 0; i < NMeasured; i++){
	for(G4int j = i+1; j < NMeasured; j++){

	  if( ( (measuredX[i] - measuredX[j])*(measuredX[i] - measuredX[j])
		   + (measuredY[i] - measuredY[j])*(measuredY[i] - measuredY[j])
		   + (measuredZ[i] - measuredZ[j])*(measuredZ[i] - measuredZ[j])
		   < packingRes2 )                                  // proximal
	      && detNum[i] == detNum[j]                          // same crystal
	      && segNum[i] == segNum[j]                          // same segment
	      &&  (NCons[i] > 0 && NCons[j] > 0) ){  // not already consolidated

	    // Energy-weighted average
	    if(measuredEdep[i] == 0 && measuredEdep[j] == 0){
	      measuredX[i]  = (measuredX[i] + measuredX[j])/2.;
	      measuredY[i]  = (measuredY[i] + measuredY[j])/2.;
	      measuredZ[i]  = (measuredZ[i] + measuredZ[j])/2.;

	      // Assign the earliest global time of the first raw hit in this IP
	      globalTime[i] = std::min(globalTime[i], globalTime[j]);
	    } else {
	      measuredX[i] = (measuredEdep[i]*measuredX[i] + measuredEdep[j]*measuredX[j])/(measuredEdep[i]+measuredEdep[j]);
	      measuredY[i] = (measuredEdep[i]*measuredY[i] + measuredEdep[j]*measuredY[j])/(measuredEdep[i]+measuredEdep[j]);
	      measuredZ[i] = (measuredEdep[i]*measuredZ[i] + measuredEdep[j]*measuredZ[j])/(measuredEdep[i]+measuredEdep[j]);
	      measuredEdep[i] += measuredEdep[j];

	      // Assign the earliest global time of the first raw hit in this IP
	      globalTime[i] = std::min(globalTime[i], globalTime[j]);
	    }
	    NCons[j] = -1;
	    NGammaHits--;

	  }

	}

      }

      // Calculate the total energy deposited in each segment.
      for(G4int i = 0; i < NMeasured; i++) // initialize
	if(NCons[i] > 0)
	  segmentEdep[i] = measuredEdep[i];

      G4bool singleDetector = true;
      for(G4int i = 0; i < NMeasured; i++){
	for(G4int j = i+1; j < NMeasured; j++){
	  if(NCons[i] > 0 && NCons[j] > 0
	     && detNum[i] == detNum[j]    // same crystal
	     && segNum[i] == segNum[j]){ // same segment
	    segmentEdep[i] += measuredEdep[j];
	    segmentEdep[j] += measuredEdep[i];
	  }
	  if( detNum[i] != detNum[j] )
	    singleDetector = false;
	}
      }

      // Identify events in which the full emitted gamma-ray energy
      // is deposited in a single crystal
      // (only evaluated for emitted multiplicity = 1 events).
      if( eventInfo->GetNEmittedGammas() == 1 ){

	G4double delta = totalEdep - eventInfo->GetEmittedGammaEnergy(0);

	// Threshold due to discrepancies in energy deposited by recoiling
	// Ge nuclei in pair-production events. (There are also tiny
	// discrepancies that can be positive or negative which may be
	// due to roundoff or some other error somewhere in geant4 tracking.
	// The upper bound of 0.0001 keV covers those.)
	G4double delta_th = 6.318E-5*eventInfo->GetEmittedGammaEnergy(0) + .074;

	if( singleDetector && delta > -delta_th && delta < 0.0001 )
	  eventInfo->SetFullEnergy(1);
	else
	  eventInfo->SetFullEnergy(0);

      }

      // Coordinate transformations
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
	    // (I-Yang's crmat file: translations in cm, hole indices +1.)
	    x -= crmat[h+1][c][0][3]*cm;
	    y -= crmat[h+1][c][1][3]*cm;
	    z -= crmat[h+1][c][2][3]*cm;
	    // ... then apply the inverse (transpose) rotation matrix.
	    measuredX[i] = x*crmat[h+1][c][0][0] 
	                 + y*crmat[h+1][c][1][0]
                         + z*crmat[h+1][c][2][0];
	    measuredY[i] = x*crmat[h+1][c][0][1] 
	                 + y*crmat[h+1][c][1][1]
                         + z*crmat[h+1][c][2][1];
	    measuredZ[i] = x*crmat[h+1][c][0][2] 
	                 + y*crmat[h+1][c][1][2]
                         + z*crmat[h+1][c][2][2];
	  }

	}

      }

      // Write S800 event to the output file
      if(fisInBeam)
	writeS800(timestamp, 
		  eventInfo->GetATA(), 
		  eventInfo->GetBTA(), 
		  eventInfo->GetDTA(), 
		  eventInfo->GetYTA());

      // Fold position resolution measuredX, measuredY, measuredZ
      // WARNING: this "smears" positions across the boundaries of
      //          the active volume, which doesn't correspond to
      //          the behavior of signal decomposition.
      if(posRes > 0){
	for(int i=0; i<NMeasured; i++) {
	  measuredX[i] += CLHEP::RandGauss::shoot(0, posRes);
	  measuredY[i] += CLHEP::RandGauss::shoot(0, posRes);
	  measuredZ[i] += CLHEP::RandGauss::shoot(0, posRes);
	}
      }

      // Write decomposed gamma event(s) to the output file
      writeDecomp(timestamp, 
		  NMeasured, 
		  detNum, segNum, NCons, 
		  measuredX, measuredY, measuredZ, 
		  measuredEdep, segmentEdep,
		  globalTime);
      
    } else {
      // Write S800 event to the output file with no detected gammas
      // if the allS800 flag is set.
      if(allS800 && fisInBeam)
	writeS800(timestamp, 
		  eventInfo->GetATA(), 
		  eventInfo->GetBTA(), 
		  eventInfo->GetDTA(), 
		  eventInfo->GetYTA());
    }
  }
  
  // Write emitted gamma information from this event to 
  // the output file.

  writeSim(timestamp, eventInfo);

}
// --------------------------------------------------
void EventAction::writeGEBHeader(GEBDATA* gd){

  G4int siz;

  if(print){
    GEBHEADER gh;
    gh.gd.type = gd->type;
    gh.gd.length = gd->length;
    gh.gd.timestamp = gd->timestamp;
    G4cout << "GEB header unsigned ints: ";
    for(int i=0;i<4;i++)
      G4cout << (hex) << gh.header[i] << "\t" 
	     << (dec) << gh.header[i] << G4endl;
  }

  siz = write(mode2file, (char *) gd, sizeof(GEBDATA));
  if(siz != sizeof(GEBDATA)){
    G4cout << "EventAction: error writing GEB header." << G4endl;
    G4cout << "             siz = " << siz << " (should be "
	   << sizeof(GEBDATA) << ")" << G4endl;
  }
}
// --------------------------------------------------
void EventAction::writeS800(long long int ts, G4double a, G4double b, 
			    G4double d, G4double y)
{
  if(mode2Out){
    //Construct GEB header for S800 event
    G4int siz;
    GEBDATA gd;
    gd.type = GEB_TYPE_S800PHYSDATA;
    gd.timestamp = ts;
    gd.length = sizeof(S800_PHYSICSDATA);

    //Construct payload for S800 event (mostly zeros)
    S800_PHYSICSDATA s800data;
    s800data.type =  0xABCD1234;
    s800data.crdc1_x = 0;
    s800data.crdc1_y = 0;
    s800data.crdc2_x = 0;
    s800data.crdc2_y = 0;
    s800data.ic_sum = 0;
    s800data.tof_xfp = 0;
    s800data.tof_obj = 0;
    s800data.rf = 0;
    s800data.trigger = 0;
    s800data.ic_de = 0;
    s800data.tof_xfpe1 = 0;
    s800data.tof_obje1 = 0;
    s800data.tof_rfe1 = 0;
    s800data.ata = (float)a;
    s800data.bta = (float)b;
    s800data.dta = (float)d;
    s800data.yta = (float)y;

    //Write GEB header for S800 event
    writeGEBHeader(&gd);

    //Write payload for S800 event
    siz = write(mode2file, (char *) &s800data, gd.length);
    if(siz != gd.length){
      G4cout << "EventAction: error writing S800 GEB payload." << G4endl;
      G4cout << "             siz = " << siz << " (should be "
	     << gd.length << ")" << G4endl;
    }

  }

  if(evOut && !outDetsOnly)
    evfile << "S    " 
	   << std::fixed << std::setprecision(4) 
	   << std::right << std::setw(12) 
	   << a << std::setw(12) << b << std::setw(12) 
	   << d << std::setw(12) << y << std::setw(12) 
	   << ts/10000
	   << G4endl;
}
// --------------------------------------------------
void EventAction::writeDecomp(long long int ts, 
			      G4int NMeasured, 
			      G4int detNum[], G4int segNum[], G4int NCons[],
			      G4double x[], G4double y[], G4double z[], 
			      G4double e[], G4double se[], G4double gt[])
{
  G4int siz;
  GEBDATA gd;
  CRYS_IPS crys_ips[100*MAX_INTPTS];
  G4double crys_gts[100*MAX_INTPTS][MAX_INTPTS];
  
  G4int Ndecomp = 0;
  G4bool Processed[100*MAX_INTPTS];
  for(G4int i = 0; i < NMeasured; i++)
    Processed[i] = false;

  for(G4int i = 0; i < NMeasured; i++){
    if( NCons[i] > 0 && Processed[i] == false ){

      crys_ips[Ndecomp].type = 0xABCD5678;
      crys_ips[Ndecomp].crystal_id = detNum[i]+4; // +4 to match measured data
      crys_ips[Ndecomp].num = 1;
      crys_ips[Ndecomp].tot_e = e[i]; //NOT += for the first one!

      crys_ips[Ndecomp].core_e[0] = 0;
      crys_ips[Ndecomp].core_e[1] = 0;
      crys_ips[Ndecomp].core_e[2] = 0;
      crys_ips[Ndecomp].core_e[3] = 0;
      crys_ips[Ndecomp].timestamp = ts;
      crys_ips[Ndecomp].t0 = 0;
      crys_ips[Ndecomp].cfd = 0;
      crys_ips[Ndecomp].chisq = 0;
      crys_ips[Ndecomp].norm_chisq = 0;
      crys_ips[Ndecomp].baseline = 0;
      crys_ips[Ndecomp].prestep = 0;
      crys_ips[Ndecomp].poststep = 0;
      crys_ips[Ndecomp].pad = 0;
      crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num-1 ].x = x[i];
      crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num-1 ].y = y[i];
      crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num-1 ].z = z[i];
      crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num-1 ].e = e[i];
      crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num-1 ].seg = segNum[i];
      crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num-1 ].seg_ener = se[i];
      crys_gts[Ndecomp][ crys_ips[Ndecomp].num-1 ] = gt[i];
      Processed[i] = true;

      for(G4int j = i+1; j < MAX_INTPTS; j++){ // Clear the interaction points
	  crys_ips[Ndecomp].ips[j].x        = 0.;
	  crys_ips[Ndecomp].ips[j].y        = 0.;
	  crys_ips[Ndecomp].ips[j].z        = 0.;
	  crys_ips[Ndecomp].ips[j].e        = 0.;
	  crys_ips[Ndecomp].ips[j].seg      = 0;
	  crys_ips[Ndecomp].ips[j].seg_ener = 0.;
      }

      // Get other interactions with this crystal
      for(G4int j = i+1; j < NMeasured; j++){ 
	if(NCons[j] > 0 && 
	   Processed[j] == false &&
	   detNum[j]+4 == crys_ips[Ndecomp].crystal_id){
	  crys_ips[Ndecomp].tot_e += e[j];
	  // Only MAX_INTPTS interaction points can be stored in a crys_ips.
	  if( crys_ips[Ndecomp].num < MAX_INTPTS ){
	    crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num ].x = x[j];
	    crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num ].y = y[j];
	    crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num ].z = z[j];
	    crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num ].e = e[j];
	    crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num ].seg = segNum[j];
	    crys_ips[Ndecomp].ips[ crys_ips[Ndecomp].num ].seg_ener = se[j];
	    crys_gts[Ndecomp][ crys_ips[Ndecomp].num ] = gt[j];
	  }
	  crys_ips[Ndecomp].num++; // Check for overflow below, and warn.
	  Processed[j] = true;
	}
      }
      Ndecomp++;
    }
  }

  for(G4int i = 0; i < Ndecomp; i++){

    if(timeSort){
      // Store unsorted interaction points.
      IP ips[MAX_INTPTS];
      G4double gts[MAX_INTPTS];
      std::vector<G4int> idx; // for sorted indices
      for(G4int j = 0; j < crys_ips[i].num; j++){
	idx.push_back(j);
	ips[j].x        = crys_ips[i].ips[j].x;
	ips[j].y        = crys_ips[i].ips[j].y;
	ips[j].z        = crys_ips[i].ips[j].z;
	ips[j].e        = crys_ips[i].ips[j].e;
	ips[j].seg      = crys_ips[i].ips[j].seg;
	ips[j].seg_ener = crys_ips[i].ips[j].seg_ener;
	gts[j]          = crys_gts[i][j];
      }
    
      // G4cout << "=========================" << G4endl;
      // G4cout << "Before:" << G4endl;
      // for(G4int j = 0; j < crys_ips[i].num; j++){
      //   G4cout << crys_gts[i][j] << ", " << crys_ips[i].ips[j].x
      // 	     << ", " << crys_ips[i].ips[j].y
      // 	     << ", " << crys_ips[i].ips[j].z
      // 	     << ", " << crys_ips[i].ips[j].e
      // 	     << ", " << crys_ips[i].ips[j].seg
      // 	     << ", " << crys_ips[i].ips[j].seg_ener
      // 	     << G4endl;
      // }
    
      // Time-sort the indices of the interaction points
      // (credit: https://stackoverflow.com/a/40183830)
      std::sort(idx.begin(), idx.end(), [&](int k,int l){return crys_gts[i][k] < crys_gts[i][l];} );

      // Use the time-sorted indices to re-order the interaction points.
      for(G4int j = 0; j < crys_ips[i].num; j++){
	crys_ips[i].ips[j].x        = ips[idx[j]].x;
	crys_ips[i].ips[j].y        = ips[idx[j]].y;
	crys_ips[i].ips[j].z        = ips[idx[j]].z;
	crys_ips[i].ips[j].e        = ips[idx[j]].e;
	crys_ips[i].ips[j].seg      = ips[idx[j]].seg;
	crys_ips[i].ips[j].seg_ener = ips[idx[j]].seg_ener;
	crys_gts[i][j]              = gts[idx[j]];
      }

      // G4cout << "After:" << G4endl;
      // for(G4int j = 0; j < crys_ips[i].num; j++)
      //   G4cout << crys_gts[i][j] << ", " << crys_ips[i].ips[j].x
      // 	     << ", " << crys_ips[i].ips[j].y
      // 	     << ", " << crys_ips[i].ips[j].z
      // 	     << ", " << crys_ips[i].ips[j].e
      // 	     << ", " << crys_ips[i].ips[j].seg
      // 	     << ", " << crys_ips[i].ips[j].seg_ener
      // 	     << G4endl;

    }
    
    if(crys_ips[i].num > MAX_INTPTS){
      G4cerr << "Warning: " << crys_ips[i].num << " interaction points."
	     << "         only " << MAX_INTPTS << " can be written."
	     << " (event " << ts/10000 << ")"
	     << G4endl;
      crys_ips[i].num = MAX_INTPTS;
    }
  }

  
  if(mode2Out){
    //Construct GEB header for decomp event(s)
    gd.type = GEB_TYPE_DECOMP;
    gd.timestamp = ts;
    gd.length = sizeof(CRYS_IPS);

    for(G4int i = 0; i < Ndecomp; i++){

      //Write GEB header for decomp event
      writeGEBHeader(&gd);

      //Write GEB payload for decomp event
      siz = write(mode2file, (char *) &crys_ips[i], sizeof(CRYS_IPS));
      if(siz != sizeof(CRYS_IPS)){
	G4cout << "EventAction: error writing decomp GEB payload." << G4endl;
	G4cout << "             siz = " << siz << " (should be "
	       << sizeof(CRYS_IPS) << ")" << G4endl;
      }

    }

  }
  
  if(evOut){
    evfile << "D" << std::setw(4) << Ndecomp 
    	   << std::setw(12) << ts/10000 << G4endl;
    for(G4int i = 0; i < Ndecomp; i++){
      evfile << "C" << std::setw(4) << crys_ips[i].crystal_id
	     << std::setw(4) << crys_ips[i].num << G4endl;
      for(G4int j = 0; j < crys_ips[i].num; j++){
	evfile << std::setw(5) 
	       << crys_ips[i].ips[j].seg 
	       << std::fixed << std::setprecision(4) 
	       << std::right << std::setw(12) 
	       << crys_ips[i].ips[j].e << std::setw(12) 
	       << crys_ips[i].ips[j].x << std::setw(12) 
	       << crys_ips[i].ips[j].y << std::setw(12) 
	       << crys_ips[i].ips[j].z << std::setw(12)
	       << crys_gts[i][j]
	       << G4endl;
      }
    }
  }

}
// --------------------------------------------------
void EventAction::writeSim(long long int ts, EventInformation* eventInfo)
{
  if(mode2Out){
    G4int siz;
    GEBDATA gd;
    G4SIM_EGS g4sim_egs;

    //Construct GEB header for G4SIM event
    gd.type = GEB_TYPE_G4SIM;
    gd.timestamp = ts;
    gd.length = sizeof(G4SIM_EGS);

    //Construct GEB payload for G4SIM event
    g4sim_egs.type = 0xABCD1234;
    g4sim_egs.num = eventInfo->GetNEmittedGammas();
    g4sim_egs.full = eventInfo->GetFullEnergy();

    for(G4int i = 0; i < g4sim_egs.num; i++){
      g4sim_egs.gammas[i].e     = eventInfo->GetEmittedGammaEnergy(i);
      g4sim_egs.gammas[i].x     = eventInfo->GetEmittedGammaPosX(i);
      g4sim_egs.gammas[i].y     = eventInfo->GetEmittedGammaPosY(i);
      g4sim_egs.gammas[i].z     = eventInfo->GetEmittedGammaPosZ(i);
      g4sim_egs.gammas[i].phi   = eventInfo->GetEmittedGammaPhi(i);
      g4sim_egs.gammas[i].theta = eventInfo->GetEmittedGammaTheta(i);
      g4sim_egs.gammas[i].beta  = eventInfo->GetBeta(i);
    }

    //Write GEB header for G4SIM event
    writeGEBHeader(&gd);

    //Write GEB payload for G4SIM event
    siz = write(mode2file, (char *) &g4sim_egs, sizeof(G4SIM_EGS));
    if(siz != sizeof(G4SIM_EGS)){
      G4cout << "EventAction: error writing G4SIM GEB payload." << G4endl;
      G4cout << "             siz = " << siz << " (should be " 
	     << sizeof(G4SIM_EGS) << ")" << G4endl;
    }

  }

  if(evOut && !outDetsOnly){
    evfile << "E" << std::setw(4) << eventInfo->GetNEmittedGammas()  
	   << std::setw(4) << eventInfo->GetFullEnergy()  
	   << std::setw(12) << ts/10000 << G4endl;
    for(G4int i = 0; i < eventInfo->GetNEmittedGammas(); i++){
      evfile << "     "
	     << std::fixed << std::setprecision(4) 
	     << std::right << std::setw(12) 
	     << eventInfo->GetEmittedGammaEnergy(i)
	     << std::setw(12) 
	     << eventInfo->GetEmittedGammaPosX(i)
	     << std::setw(12) 
	     << eventInfo->GetEmittedGammaPosY(i)
	     << std::setw(12) 
	     << eventInfo->GetEmittedGammaPosZ(i)
	     << std::setw(12) 
	     << eventInfo->GetEmittedGammaPhi(i)
	     << std::setw(12) 
	     << eventInfo->GetEmittedGammaTheta(i)
	     << std::setw(12)
	     << eventInfo->GetBeta(i) << G4endl;
    }
  }
 
}
// --------------------------------------------------TB
void EventAction::openEvfile()
{
  if (!evfile.is_open()) evfile.open(outFileName.c_str());
  if (!evfile.is_open()){
    G4cout<< "ERROR opening evfile." << G4endl;
    evOut = false;
  } else {
    G4cout << "\nOpened output file: " << outFileName << G4endl;
    evOut = true;
  }
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
  closeEvfile();
  openEvfile();
  return;
}
// --------------------------------------------------
void EventAction::openMode2file()
{
  mode2file = open(mode2FileName.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 
		   S_IRUSR | S_IWUSR | S_IRGRP |  S_IWGRP | S_IROTH);
  if(mode2file < 0){
    G4cout << "ERROR opening mode2file" << G4endl;
    mode2Out = false;
  } else {
    G4cout << "\nOpened mode 2 Output file: " << mode2FileName << G4endl;
    mode2Out = true;
  }
  return;
}
//---------------------------------------------------
void EventAction::closeMode2file()
{
  close(mode2file);
  return;
}
//----------------------------------------------------
void EventAction::SetMode2File(G4String name)
{
  mode2FileName = name;
  openMode2file();
  return;
}
// --------------------------------------------------
void EventAction::openCrmatFile()
{
  crmatFile = open(crmatFileName.c_str(), O_RDONLY, 0);
  if(crmatFile < 0){
    G4cout << "ERROR opening crmatFile" << G4endl;
    exit(EXIT_FAILURE);
  }
  return;
}
//---------------------------------------------------
void EventAction::closeCrmatFile()
{
  close(crmatFile);
  return;
}
//---------------------------------------------------
void EventAction::SetCrmatFile(G4String name) { 

  crmatFileName = name;

  openCrmatFile();

  G4cout << "\nUsing crmat from file " << crmatFileName 
	 << " to transform interaction points from world to crytsal frames."
	 << G4endl;

  int size;
  size = read(crmatFile, (char *) crmat, sizeof(crmat));

  closeCrmatFile();

  if(size < 0){
    G4cout << "ERROR reading crmat file." << G4endl;
    exit(EXIT_FAILURE);
  }

  G4cout << "Read " << size << " bytes into crmat" << G4endl;

  if(print){
    for(int i=0;i<MAXDETPOS;i++){
      for(int j=0;j<MAXCRYSTALNO;j++){
	G4cout << "Hole : " << i << "\tCrystal: " << j << endl;
	for(int k=0;k<4;k++){	
	  for(int l=0;l<4;l++){
	    G4cout << crmat[i][j][k][l] << "\t";
	  }
	  G4cout << endl;
	} 
	G4cout <<"-------------------------------------"<< endl;
      }
      G4cout <<"-------------------------------------"<< endl;
    }
  }
  return;
}
//---------------------------------------------------
void EventAction::SetGretinaCoords(){

  gretinaCoords = true; 

  G4cout << "Writing interaction points in the Gretina coordinate system (x = down, z = beam)." << G4endl;

  return;
}
//---------------------------------------------------
 
void EventAction::SetCrystalXforms(){ 

  crystalXforms = true; 

  G4cout << "Using internal transformations from the world frame to the crystal frames for Mode 2 output." << G4endl;

  return;
}
//---------------------------------------------------
