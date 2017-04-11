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
  // G4cout << std::fixed << std::setprecision(4) << std::setw(12);
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
  if( ( aTrack->GetParticleDefinition()->GetParticleName() == "gamma" ||
	aTrack->GetParticleDefinition()->GetParticleName() == "mu+" ||
	aTrack->GetParticleDefinition()->GetParticleName() == "mu-" )
      &&
      ( aTrack->GetParentID() == 0 ||
	aTrack->GetCreatorProcess()->GetProcessName() == "Decay" ||
	aTrack->GetCreatorProcess()->GetProcessName() == "Reaction" ) ){

    //    G4cout << "> Event ID = " << eventAction->GetEvent()->GetEventID() << G4endl;

    G4ThreeVector pos = aTrack->GetPosition();
    G4ThreeVector dir = aTrack->GetMomentumDirection();

    // G4cout << std::fixed << std::setprecision(4) 
    // 	   << std::setw(12) << std::right
    // 	   << "   pos = " << pos
    // 	   << "   dir = " << dir
    // 	   << "   energy = " << aTrack->GetKineticEnergy();

    eventInfo->AddEmittedGamma(aTrack->GetKineticEnergy(), 
			       &pos, &dir,
			       aTrack->GetParentID());

  }

}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{

  eventInfo 
    = (EventInformation*)eventAction->GetEvent()->GetUserInformation();

  // G4cout << std::fixed << std::setprecision(4) << std::setw(12)
  // 	 << "> PostUserTrackingAction" << G4endl;
  // G4cout << "> KineticEnergy = " << aTrack->GetKineticEnergy() << G4endl;
  // G4cout << "> TotalEnergy = "   << aTrack->GetTotalEnergy() << G4endl;
  // G4cout << "> Pre-step Position = " 
  // 	 << aTrack->GetStep()->GetPreStepPoint()->GetPosition() << G4endl;
  // G4cout << "> Post-step Position = " 
  // 	 << aTrack->GetStep()->GetPostStepPoint()->GetPosition() << G4endl;
  // G4cout << "> MomentumDirection = " 
  // 	 << aTrack->GetMomentumDirection() << G4endl;
  // G4cout << "> Beta = " 
  // 	 << aTrack->GetStep()->GetPostStepPoint()->GetBeta() << G4endl;
  // G4cout << ">==========================================" << G4endl;

  // S800 data
  //
  // Geant4 coordinates:   x: North, y: up,    z: with the beam
  //
  // NSCL coordinates:     x: down,  y: North, z: with the beam
  //                  (x: ~up in the S800 focal plane detectors)

  if( aTrack->GetParticleDefinition()->GetParticleType() == "nucleus" &&
      aTrack->GetParentID() > 0 ){

    if( !aTrack->GetDefinition ()-> GetPDGStable() )
      eventInfo->AddBeta(aTrack->GetStep()->GetPostStepPoint()->GetBeta(),
			 aTrack->GetTrackID());
    
  }

}
