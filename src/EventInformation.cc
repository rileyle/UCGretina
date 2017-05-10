#include "EventInformation.hh"

EventInformation::EventInformation() { 
  fNEmittedGammas = 0;
  fNBetas = 0;
  fFullEnergy = -1;
  fata      = sqrt(-1.0);
  fbta      = sqrt(-1.0);
  fdta      = sqrt(-1.0);
  fyta      = sqrt(-1.0);
}

void EventInformation::AddEmittedGamma(G4double e, 
				       G4ThreeVector *pos, 
				       G4ThreeVector *dir,
				       G4int parentID){

  //G4cout << "   fNEmittedGammas = " << fNEmittedGammas << G4endl;

  for(G4int i = 0; i < fNBetas; i++){
    // G4cout << "i have fBetaTrackID[" <<i<<"] = "<<fBetaTrackID[i]
    // 	   << ", and Emitted Gamma parent ID = " << parentID << G4endl;
    if( fBetaTrackID[i] == parentID ){
      // G4cout << "... saving energy "<<e/keV<<" keV with beta = "
      // 	     <<fBeta[i]<<G4endl;
      fEmittedGammaEnergies[i] = e/keV;
      fEmittedGammaPosX[i]     = pos->getX()/mm;
      fEmittedGammaPosY[i]     = pos->getY()/mm;
      fEmittedGammaPosZ[i]     = pos->getZ()/mm;
      fEmittedGammaPhi[i]      = dir->getPhi()/rad;
      fEmittedGammaTheta[i]    = dir->getTheta()/rad;
      fNEmittedGammas++;
      break;
    }
  }
}

void EventInformation::AddBeta(G4double b, G4int TID){

  // Beta values are stored first, so we keep track of the track ID 
  // of the ion to match it with the corresponding emitted gamma
  // in AddEmittedGamma.
  
  fBeta[fNBetas] = b;
  fBetaTrackID[fNBetas] = TID;
  fNBetas++;

}
