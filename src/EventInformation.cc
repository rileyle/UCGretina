#include "EventInformation.hh"

EventInformation::EventInformation() { 
  fNEmittedGammas = 0;
  fata      = sqrt(-1.0);
  fbta      = sqrt(-1.0);
  fdta      = sqrt(-1.0);
  fyta      = sqrt(-1.0);
}

void EventInformation::AddEmittedGamma(G4double e, 
				       G4ThreeVector *pos, 
				       G4ThreeVector *dir){

  //  G4cout << "   fNEmittedGammas = " << fNEmittedGammas << G4endl;
  fEmittedGammaEnergies[fNEmittedGammas] = e/keV;
  fEmittedGammaPosX[fNEmittedGammas] = pos->getX()/mm;
  fEmittedGammaPosY[fNEmittedGammas] = pos->getY()/mm;
  fEmittedGammaPosZ[fNEmittedGammas] = pos->getZ()/mm;
  fEmittedGammaPhi[fNEmittedGammas] = dir->getPhi()/rad;
  fEmittedGammaTheta[fNEmittedGammas] = dir->getTheta()/rad;
  fNEmittedGammas++;

}
