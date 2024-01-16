//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "EventAction.hh"
#include "EventInformation.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction()                                         
{
  eventAction = (EventAction*)
    G4RunManager::GetRunManager()->GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  eventInfo 
    = (EventInformation*)eventAction->GetEvent()->GetUserInformation();

  // if(   aStep->GetTrack()->GetDefinition()->GetParticleType() == "nucleus" 
  // 	&& aStep->GetPostStepPoint()->GetStepStatus() != fWorldBoundary ){
  // 	//  	&& aStep->GetTrack()->GetCurrentStepNumber() > 30){

  //   G4cout << std::fixed << std::setprecision(4)
  //   	   << "event: " << eventAction->GetEvent()->GetEventID()
  // 	   << ", " << aStep->GetTrack()->GetDefinition()->GetParticleName()
  //   	   << ", parentID: " << aStep->GetTrack()->GetParentID()
  //   	   << ", v1: " << aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
  //   	   << ", v2: " << aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
  //   	   << ", KE: " << aStep->GetTrack()->GetKineticEnergy() 
  //   	   << ", pos1 (" 
  //   	   << aStep->GetPreStepPoint()->GetPosition().getX()
  //   	   << ", "
  //   	   << aStep->GetPreStepPoint()->GetPosition().getY()
  //   	   << ", "
  //   	   << aStep->GetPreStepPoint()->GetPosition().getZ()
  //   	   << ")" 
  //   	   << ", pos2 (" 
  //   	   << aStep->GetPostStepPoint()->GetPosition().getX()
  //   	   << ", "
  //   	   << aStep->GetPostStepPoint()->GetPosition().getY()
  //   	   << ", "
  //   	   << aStep->GetPostStepPoint()->GetPosition().getZ()
  //   	   << ")" 
  // 	   << " step " << aStep->GetTrack()->GetCurrentStepNumber()
  // 	   << G4endl;
  // }

  G4String particleName = aStep->GetTrack()->GetDefinition()->GetParticleName();
  
  if( aStep->GetTrack()->GetDefinition()->GetParticleType() == "nucleus" 
      && ( aStep->GetTrack()->GetParentID() > 0 || eventAction->AllS800() )
      && aStep->GetPostStepPoint()->GetStepStatus() != fWorldBoundary ){

    // get initial and final volumes of the current step
    G4VPhysicalVolume* volume1 
      = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
    G4VPhysicalVolume* volume2
      = aStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume();

    // Record S800 data as the reaction product leaves the target
    if(   G4StrUtil::contains(volume1->GetName(), "Target")
       && G4StrUtil::contains(volume2->GetName(), "expHall") ){

      G4ThreeVector pDir = aStep->GetTrack()->GetMomentumDirection();

      // ATA is the dispersive angle (radians), 
      // down is + in NSCL coords= -y in Geant4 coords
      eventInfo->SetATA( asin(-pDir.getY()/pDir.mag()) );

      // BTA is the non-dispersive angle (radians), 
      // South is + in NSL coords = -x in Geant4 coords
      eventInfo->SetBTA( asin(-pDir.getX()/pDir.mag()) );

      // DTA is dT/T with T = kinetic energy corresponding 
      // to the user-supplied center of the S800 acceptance
      eventInfo->SetDTA( (aStep->GetTrack()->GetKineticEnergy() 
			  - eventAction->GetS800KE())
			  / eventAction->GetS800KE() ); 

      // YTA is horizontal position on target (mm), 
      // South is + in NSCL coords = -x in Geant4 coords
      eventInfo->SetYTA( -aStep->GetTrack()->GetStep()->GetPreStepPoint()->GetPosition().getX()*mm );

    }
    // Kill a reaction product once it hits the chamber or beamline
    // as long it has already emitted its gamma(s)
    else if( ( G4StrUtil::contains(volume1->GetName(), "BeamTube")
	       || G4StrUtil::contains(volume1->GetName(), "BeamTee")
	       || G4StrUtil::contains(volume1->GetName(), "Cell")
	       || G4StrUtil::contains(volume1->GetName(), "Bar")
	       || G4StrUtil::contains(volume1->GetName(), "Chamber")
	       || G4StrUtil::contains(volume1->GetName(), "sled")
	       || G4StrUtil::contains(volume1->GetName(), "Frame")
	       || G4StrUtil::contains(volume1->GetName(), "Tape")
	       || G4StrUtil::contains(volume1->GetName(), "Ring") )
	     && !G4StrUtil::contains(aStep->GetTrack()->GetParticleDefinition()->GetParticleName(), '[') ){

           // G4cout << "************************* SteppingAction: terminating track in "
      	   //   << volume1->GetName()
      	   //   << endl;

      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
	
    }
    // Kill light charged particles at rest
    else if( ( particleName == "proton"    || particleName == "deuteron"
	       || particleName == "triton" || particleName == "alpha" )
	     && aStep->GetTrack()->GetKineticEnergy() < 0.1*keV ){

      aStep->GetTrack()->SetTrackStatus(fStopAndKill);
      
    }

  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
