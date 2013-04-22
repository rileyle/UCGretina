#include "TrackingAction.hh"

TrackingAction::TrackingAction(EventAction* evt){
  eventAction = evt;
}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // G4cout << "\n>==========================================" << G4endl;
  // G4cout << "> PreUserTrackingAction" << G4endl;
  // G4cout << "> Event ID = " << eventAction->GetEvent()->GetEventID() << G4endl;
  // G4cout << "> TrackID = "       << aTrack->GetTrackID() << G4endl;
  // G4cout << "> KineticEnergy = " << aTrack->GetKineticEnergy() << G4endl;
  // G4cout << "> TotalEnergy = "   << aTrack->GetTotalEnergy() << G4endl;
  // G4cout << "> Position = "      << aTrack->GetPosition() << G4endl;
  // G4cout << "> MomentumDirection = " 
  // 	 << aTrack->GetMomentumDirection() << G4endl;
  // G4cout << "> ParticleDefinition Name = "
  // 	 << aTrack->GetParticleDefinition()->GetParticleName() << G4endl;
  // G4cout << "> ParticleDefinition Type = "
  // 	 << aTrack->GetParticleDefinition()->GetParticleType() << G4endl;
  // G4cout << "> Parent ID = " 
  // 	 << aTrack->GetParentID() << G4endl;
  // if(aTrack->GetParentID() > 0)
  //   G4cout << "> CreatorProcess Name = " 
  // 	   << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
  // G4cout << ">      ----------------------------" << G4endl;

  eventInfo 
    = (EventInformation*)eventAction->GetEvent()->GetUserInformation();

  // Emitted gamma
  if( aTrack->GetParticleDefinition()->GetParticleName() == "gamma" &&
      ( aTrack->GetParentID() ==0 ||
	aTrack->GetCreatorProcess()->GetProcessName() == "Decay" ) ){
    G4ThreeVector pos = aTrack->GetPosition();
    G4ThreeVector dir = aTrack->GetMomentumDirection();
    eventInfo->AddEmittedGamma(aTrack->GetKineticEnergy(), 
			       &pos, &dir);
  }

}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  // G4cout << "> PostUserTrackingAction" << G4endl;
  // G4cout << "> KineticEnergy = " << aTrack->GetKineticEnergy() << G4endl;
  // G4cout << "> TotalEnergy = "   << aTrack->GetTotalEnergy() << G4endl;
  // G4cout << "> Pre-step Position = " 
  // 	 << aTrack->GetStep()->GetPreStepPoint()->GetPosition() << G4endl;
  // G4cout << "> Post-step Position = " 
  // 	 << aTrack->GetStep()->GetPostStepPoint()->GetPosition() << G4endl;
  // G4cout << "> MomentumDirection = " 
  // 	 << aTrack->GetMomentumDirection() << G4endl;
  // G4cout << ">==========================================" << G4endl;

  // S800 data
  //
  // Geant4 coordinates:   x: North, y: up,    z: with the beam
  //
  // NSCL coordinates:     x: down,  y: North, z: with the beam
  //                  (x: ~up in the S800 focal plane detectors)

  if( aTrack->GetParticleDefinition()->GetParticleType() == "nucleus" &&
      aTrack->GetParentID() > 0 ){
    if( aTrack->GetCreatorProcess()->GetProcessName() == "Decay" ) {
      G4ThreeVector pDir = aTrack->GetMomentumDirection();
      // ATA is the dispersive angle, down is + in NSCL coords= -y in Geant4 coords
      eventInfo->SetATA( asin(-pDir.getY()/pDir.mag())/mrad );
      // BTA is the non-dispersive angle, South is + in NSL coords = -x in Geant4 coords
      eventInfo->SetBTA( asin(-pDir.getX()/pDir.mag())/mrad );
      // DTA is dT/T with T = kinetic energy corresponding to the center of the S800 acceptance
      eventInfo->SetDTA((aTrack->GetKineticEnergy() - eventAction->GetS800KE()) / eventAction->GetS800KE() ); 
    } else if ( aTrack->GetCreatorProcess()->GetProcessName() == "Reaction" ) {
      // YTA is horizontal position on target, South is + in NSCL coords = -x in Geant4 coords
      eventInfo->SetYTA( -aTrack->GetStep()->GetPreStepPoint()->GetPosition().getX()/mm );
    }
  }

}
