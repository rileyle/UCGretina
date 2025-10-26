#include "PrimaryVertexInformation.hh"

PrimaryVertexInformation::PrimaryVertexInformation() { 
  fNEmittedGammas = 0;
  fNBetas         = 0;
  fFullEnergy     = -1;
  fata            = sqrt(-1.0);
  fbta            = sqrt(-1.0);
  fdta            = sqrt(-1.0);
  fyta            = sqrt(-1.0);
  fExitPos        = G4ThreeVector(0, 0, 0);
  fExitBeta       = sqrt(-1.0);
  fExitTheta      = sqrt(-1.0);
  fExitPhi        = sqrt(-1.0);
  ffiltercode     = 0;
  fwrite          = true;
}

void PrimaryVertexInformation::AddEmittedGamma(G4double e, 
					       G4ThreeVector *pos, 
					       G4ThreeVector *dir,
					       G4int parentID){

  //G4cout << "   fNEmittedGammas = " << fNEmittedGammas << G4endl;

  // In-beam 
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

  // Source/Background/Cache
  if(fNBetas == 0){
      fEmittedGammaEnergies[fNEmittedGammas] = e/keV;
      fEmittedGammaPosX[fNEmittedGammas]     = pos->getX()/mm;
      fEmittedGammaPosY[fNEmittedGammas]     = pos->getY()/mm;
      fEmittedGammaPosZ[fNEmittedGammas]     = pos->getZ()/mm;
      fEmittedGammaPhi[fNEmittedGammas]      = dir->getPhi()/rad;
      fEmittedGammaTheta[fNEmittedGammas]    = dir->getTheta()/rad;
      fBeta[fNBetas] = 0;
      fNEmittedGammas++;
  }
    
}

void PrimaryVertexInformation::AddBeta(G4double b, G4int TID){

  // Beta values are stored first, so we keep track of the track ID 
  // of the ion to match it with the corresponding emitted gamma
  // in AddEmittedGamma.
  
  fBeta[fNBetas] = b;
  fBetaTrackID[fNBetas] = TID;
  fNBetas++;

}

void PrimaryVertexInformation::SetExitPos(G4ThreeVector* p){
  fExitPos.setX(p->x());
  fExitPos.setY(p->y());
  fExitPos.setZ(p->z());
}
