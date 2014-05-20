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
				       G4ThreeVector *dir){

  //G4cout << "   fNEmittedGammas = " << fNEmittedGammas << G4endl;
  fEmittedGammaEnergies[fNEmittedGammas] = e/keV;
  fEmittedGammaPosX[fNEmittedGammas] = pos->getX()/mm;
  fEmittedGammaPosY[fNEmittedGammas] = pos->getY()/mm;
  fEmittedGammaPosZ[fNEmittedGammas] = pos->getZ()/mm;
  fEmittedGammaPhi[fNEmittedGammas] = dir->getPhi()/rad;
  fEmittedGammaTheta[fNEmittedGammas] = dir->getTheta()/rad;
  fNEmittedGammas++;

}

void EventInformation::AddBeta(G4double b){

  fBeta[fNBetas] = b;
  fNBetas++;

}
