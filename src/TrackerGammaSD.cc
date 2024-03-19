
#include "TrackerGammaSD.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4VTouchable.hh"
#include <string.h>
#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerGammaSD::TrackerGammaSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="gammaCollection");
  print = false;
  phdA = 0.21;
  phdB = 1.099;
  posRes = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerGammaSD::~TrackerGammaSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerGammaSD::Initialize(G4HCofThisEvent*)
{

  gammaCollection = new TrackerGammaHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 

  G4RunManager* runManager = G4RunManager::GetRunManager();
  if(runManager->GetCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetG4code()->GetParticleName()
     == "gamma")
    primaryGammaID = 0; // Source and background simulations
  else
    primaryGammaID = 2; // In-beam simulations

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerGammaSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();

  G4String name = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  G4int detCode = -1;
  G4int detNum = -1;
  G4VSolid* solid = theTouchable->GetSolid(0);
  if(name == "LaBr"){
    detCode = 0;
    detNum  = 132;
  }
  else if(name.contains("Leaf")){ // Clover
    detCode = theTouchable->GetReplicaNumber(0);
    detNum = detCode%1000;
  }
  else { // GRETINA
    detCode = theTouchable->GetReplicaNumber(0);
    if(detCode == 0){
      detCode = theTouchable->GetReplicaNumber(depth);
      solid = theTouchable->GetSolid(depth);
    }
    detNum = detCode%1000;
  }

  //  G4cout << name << ", " << detCode << ", " << detNum << G4endl;

  // Track everything (the primary gamma, all secondary electrons, 
  // and all secondary gammas).
  G4double edep = aStep->GetTotalEnergyDeposit();

  G4String particleType = aStep->GetTrack()->GetDefinition()->GetParticleType();
  G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
  G4String processName
    = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  //  G4cout << " &&& particleType = " << particleType 
  //	 << ", particleName = " << particleName << G4endl;
  
  // Pulse height defect for (n,n')
  // Joa Ljungvall and Johan Nyberg, NPA546, 553–573 (2005), Figure 3
  if( ( particleType == "nucleus" && processName == "ionIoni" ) ||
      ( particleName == "neutron" && processName == "hadElastic" ) )
    edep = phdA*pow(edep, phdB);

  // Keep the gamma parent of pair-production and polarized Compton tracks for hit processing.
  if(edep<0.001*eV
     && processName != "conv" && processName != "pol-compt" && processName != "pol-conv") 
    return false;

  G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
  
  G4VPhysicalVolume* topVolume;
  if( detNum < 124 ) // GRETINA
    topVolume = theTouchable->GetVolume(depth);
  else if( detNum > 123 ) // Clover or LaBr
    topVolume = theTouchable->GetVolume(0);
  else
    topVolume = NULL;

  G4ThreeVector frameTrans = topVolume->GetFrameTranslation();

  const G4RotationMatrix* rot = topVolume->GetFrameRotation();
                                  
  G4RotationMatrix frameRot;
  if( rot )
    frameRot = *rot;
  
  G4ThreeVector posSol = frameRot( position );
  posSol += frameRot( frameTrans );

  // "Measured" position, with simulated position resolution
  // WARNING: this "smears" positions across the boundaries of
  //          the active volume, which doesn't correspond to
  //          the behavior of signal decomposition.
  G4ThreeVector positionM = position;
  if(posRes > 0){
    positionM.setX(position.getX() + CLHEP::RandGauss::shoot(0, posRes));
    positionM.setY(position.getY() + CLHEP::RandGauss::shoot(0, posRes));
    positionM.setZ(position.getZ() + CLHEP::RandGauss::shoot(0, posRes));
    G4ThreeVector posSolM = frameRot( positionM );
    posSolM += frameRot( frameTrans );
    if(solid->Inside(posSolM) == kOutside){
      // Displacement in Solid (Crystal-frame) coordinates
      G4ThreeVector deltaSol = posSolM-posSol;
      // Distance to the surface along the displacement vector.
      G4double dr = solid->DistanceToOut(posSol, deltaSol/deltaSol.mag());
      // Displacement in world coordinates
      G4ThreeVector delta = positionM - position;
      // Corrected "smeared" position (at surface)
      positionM = position + delta/delta.mag()*dr;
      //      posSolM = frameRot( positionM );
      //      posSolM += frameRot( frameTrans );
      //      G4cout <<"new insideM = " << solid->Inside(posSolM) << G4endl;
    }
  } else
    positionM.setX(1.0e12); // Block processing/printing of smeared positions 
                            // if the position resolution is set to 0.
  
  // Segment number
  G4int segCode = 0;
  G4RunManager* runManager = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector 
    = (DetectorConstruction*)runManager->GetUserDetectorConstruction();

  if(theDetector->GetGretina()->GetReadOut()){
    if( detNum < 124 ){ // GRETINA
      segCode = 
	theDetector->GetGretina()->GetSegmentNumber( detCode, posSol );

      // Modify sector number to match GRETINA data stream
      G4int slice  = segCode/10;
      G4int sector = segCode%10;
      // Type B crystal (offset -1)
      if(sector>0) 
	sector--;
      else
	sector=5;
      // Type A crystal (offset -2)
      if(detNum % 2){
	if(sector>0) 
	  sector--;
	else
	  sector=5;
      }
      segCode = sector + 6 * slice;
    }
  } else if( name.contains("Leaf") ){ // Clover
    // Crystal/segment labels follow
    // Eurysys CLOVER 4X50X80 SEG2 manual p. 18
    if( name.contains("pv_0") ){ // Crystal 1
      if(posSol.getX() > 0)
    	segCode = 1; // L
      else if(posSol.getX() < 0)
    	segCode = 2; // M
    }
    if( name.contains("pv_1") ){ // Crystal 2
      if(posSol.getY() > 0)
    	segCode = 3; // R
      else if(posSol.getY() < 0)
    	segCode = 2; // M
    }
    if( name.contains("pv_2") ){ // Crystal 3
      if(posSol.getX() > 0)
    	segCode = 3; // R
      else if(posSol.getX() < 0)
    	segCode = 2; // M
    }
    if( name.contains("pv_3") ){ // Crystal 4
      if(posSol.getY() > 0)
    	segCode = 1; // L
      else if(posSol.getY() < 0)
    	segCode = 2; // M
    }
  } else {
    segCode = -1;
  }
  TrackerGammaHit* newHit = new TrackerGammaHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());

  newHit->SetParticleID(particleName);
  newHit->SetProcess(processName);
  newHit->SetParentTrackID(aStep->GetTrack()->GetParentID());
  if( aStep->GetTrack()->GetCreatorProcess() )
    newHit->SetCreatorProcess(aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
  else
    newHit->SetCreatorProcess("primary");
  newHit->SetDetNumb    (detNum);
  newHit->SetSegNumb    (segCode);
  newHit->SetEdep       (edep);
  newHit->SetKE         (aStep->GetTrack()->GetKineticEnergy());
  newHit->SetPos        (position);
  newHit->SetPosM       (positionM);
  newHit->SetPosCrys    (posSol);
  newHit->SetTrackOrigin(aStep->GetTrack()->GetVertexPosition());
  newHit->SetGlobalTime(aStep->GetPostStepPoint()->GetGlobalTime());
  
  gammaCollection->insert( newHit );
  newHit->Draw();
  //  getc(stdin);
  //  newHit->Print(); 
  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerGammaSD::EndOfEvent(G4HCofThisEvent* HCE)
{
   G4int i;
   G4int NbHits = gammaCollection->entries();

       if (NbHits>0&&print)

	  { 

	    G4RunManager* runManager = G4RunManager::GetRunManager();
	    
	    G4cout << "\n--------> event " << runManager->GetCurrentEvent()->GetEventID() << ", "
		   << NbHits << " hits for gamma tracking: " << G4endl;
	    G4cout << "                            parent    creator" << G4endl;
	    G4cout << "trackID   PID     process   track     process      det seg     Edep      KE         X         Y         Z         Xo        Yo        Zo    T (ps)" << G4endl;
	    G4double totE = 0;
	    for (i=0;i<NbHits;i++){
	      (*gammaCollection)[i]->Print();
	      totE += (*gammaCollection)[i]->GetEdep();
	    }
	    G4cout << "total energy deposited = "
		   << std::fixed << std::setprecision(2)
		   << std::setw(10) << std::right
		   << totE/keV << " keV" << G4endl;
	    
	  }


  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, gammaCollection ); 
 }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

