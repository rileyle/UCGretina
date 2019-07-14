
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
  if(name == "LaBr"){
    detCode = 0;
    detNum  = 132;
  }
  else { // GRETINA or Clover
    detCode = theTouchable->GetReplicaNumber(0);
    if(detCode == 0)
      detCode = theTouchable->GetReplicaNumber(depth);
    detNum = detCode%1000;
  }

  //  G4cout << name << ", " << detCode << ", " << detNum << G4endl;

  // Track everything (the primary gamma, all secondary electrons, 
  // and all secondary gammas).
  G4double edep = aStep->GetTotalEnergyDeposit();

  // Keep the gamma parent of pair-production tracks for hit processing.
  if(edep<0.001*eV
     && aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() 
     != "conv") 
    return false;

  G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();

  G4VPhysicalVolume* topVolume;
  //#ifdef SCANNING
  if( detNum < 124 ) // GRETINA
  //#endif
    topVolume = theTouchable->GetVolume(depth);
  //#ifdef SCANNING
  else if( detNum > 123 ) // Clover or LaBr
    topVolume = theTouchable->GetVolume(0);
  else
    topVolume = NULL;
  //#endif

  G4ThreeVector frameTrans = topVolume->GetFrameTranslation();

  const G4RotationMatrix* rot = topVolume->GetFrameRotation();
                                  
  G4RotationMatrix frameRot;
  if( rot )
    frameRot = *rot;
  
  G4ThreeVector posSol = frameRot( position );
  posSol += frameRot( frameTrans );

  // Segment number
  G4int segCode = 0;
  G4RunManager* runManager = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector 
    = (DetectorConstruction*)runManager->GetUserDetectorConstruction();

  if(theDetector->GetGretina()->GetReadOut()){
    //#ifdef SCANNING
    if( detNum < 124 ){ // GRETINA
    //#endif
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
    //#ifdef SCANNING
    }
    //#endif
  } else {
    segCode = -1;
  }
  TrackerGammaHit* newHit = new TrackerGammaHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());

  newHit->SetParticleID(aStep->GetTrack()->GetDefinition()->GetParticleName());
  newHit->SetProcess(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
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
  newHit->SetPosCrys    (posSol);
  newHit->SetTrackOrigin(aStep->GetTrack()->GetVertexPosition());

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
	    G4cout << "trackID   PID     process   track     process     det seg     Edep      KE         X         Y         Z         Xo        Yo        Zo" << G4endl;
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

