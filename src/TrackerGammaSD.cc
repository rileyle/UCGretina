
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

  G4String name = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if(name.substr(0,5)=="geCap")
    {

      G4RunManager* runManager = G4RunManager::GetRunManager();
      DetectorConstruction* theDetector = (DetectorConstruction*)runManager->GetUserDetectorConstruction();

      // Track everything (the primary gamma, all secondary electrons, and all secondary gammas).
      G4double edep = aStep->GetTotalEnergyDeposit();

      // Keep the gamma parent of pair-production tracks for hit processing.
      if(edep<0.001*eV
	 && aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() != "conv") 
	return false;

      G4int detCode = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(depth);
      G4int detNum  = detCode%1000;

      G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();

      ////////////////////////////////////////////////////////////////////
      /// Position of the interaction point in the solid reference frame
      ///////////////////////////////////////////////////////////////////
      G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
      G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
      G4VPhysicalVolume* topVolume = theTouchable->GetVolume(depth);
  
      G4ThreeVector frameTrans = topVolume->GetFrameTranslation();
                                  
      const G4RotationMatrix* rot = topVolume->GetFrameRotation();
                                  
      G4RotationMatrix frameRot;
      if( rot )
	frameRot = *rot;
  
      G4ThreeVector posSol = frameRot( position );
      posSol += frameRot( frameTrans );
    
      //////////////////////
      /// Segment number
      //////////////////////
      G4int segCode = 0;

      segCode = theDetector->GetGretina()->GetSegmentNumber( 0, detCode, posSol );

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
      newHit->SetPos        (position);
      newHit->SetTrackOrigin(aStep->GetTrack()->GetVertexPosition());

      gammaCollection->insert( newHit );
      newHit->Draw();
      //  getc(stdin);
  //   newHit->Print(); 
      return true;
      }
  else
    {
      //    G4cout << "Energy deposit in the " << name << G4endl;
      //    G4cout << "E="<<G4BestUnit(edep,"Energy")<<G4endl;
      //    G4cout <<"Event ignored" <<G4endl;
	    //	    getc(stdin);
      return false;
    }

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

	    for (i=0;i<NbHits;i++) (*gammaCollection)[i]->Print();

	  }


  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, gammaCollection ); 
 }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

