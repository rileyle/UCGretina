
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
  print=false; //LR (formerly not initialized)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerGammaSD::~TrackerGammaSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerGammaSD::Initialize(G4HCofThisEvent*)
{


  gammaCollection = new TrackerGammaHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerGammaSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4String name = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if(name.substr(0,5)=="geCap")
    {
      G4RunManager* runManager = G4RunManager::GetRunManager();
      DetectorConstruction* theDetector = (DetectorConstruction*)runManager->GetUserDetectorConstruction();

      G4double DE = aStep->GetTotalEnergyDeposit();
      if(DE<0.001*eV) return false;

      G4double edep = aStep->GetTotalEnergyDeposit();

      G4double etotal = aStep->GetTrack()->GetVertexKineticEnergy(); 

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

      TrackerGammaHit* newHit = new TrackerGammaHit();

      newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());

      newHit->SetParticleID(aStep->GetTrack()->GetDefinition()->GetParticleName());
      newHit->SetProcess(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());

      newHit->SetDetNumb(detNum);
      newHit->SetSegNumb(segCode);
      newHit->SetEdep     (edep);
      newHit->SetTotalEnergy(etotal);
      newHit->SetPos      (position);

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
	    
	   G4cout << "\n--------> " << NbHits << " hits for gamma tracking: " << G4endl;

	   for (i=0;i<NbHits;i++) (*gammaCollection)[i]->Print();

	   }


  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, gammaCollection ); 
 }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

