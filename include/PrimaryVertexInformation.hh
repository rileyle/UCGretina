#ifndef PrimaryVertexInformation_h
#define PrimaryVertexInformation_h 1

#include "G4VUserPrimaryVertexInformation.hh"
#include "G4ThreeVector.hh"
#include "GEB.hh"
#include "G4SystemOfUnits.hh"

class PrimaryVertexInformation : public G4VUserPrimaryVertexInformation
{
public:

  PrimaryVertexInformation();
  ~PrimaryVertexInformation(){};

  inline virtual void Print()const{;}

  void AddEmittedGamma(G4double, G4ThreeVector*, G4ThreeVector*, G4int);
  void AddBeta(G4double, G4int);
  void SetFullEnergy(G4int f){ fFullEnergy = f; }
  void SetATA(G4double a){fata = a;}
  void SetBTA(G4double b){fbta = b;}
  void SetDTA(G4double d){fdta = d;}
  void SetYTA(G4double y){fyta = y;}
  void SetExitPos(G4ThreeVector* p);
  void SetExitBeta(G4double b){ fExitBeta = b; }
  void SetExitTheta(G4double t){ fExitTheta = t; }
  void SetExitPhi(G4double p){ fExitPhi = p; }
  void SetExitTime(G4double t){ fExitTime = t; }
  void SetFilterCode(G4int c){ffiltercode = c;}
  void SetWriteEvent(G4bool b){fwrite = b;}

  G4double GetEmittedGammaEnergy(G4int i){ return fEmittedGammaEnergies[i]; }
  G4double GetEmittedGammaPosX(G4int i){ return fEmittedGammaPosX[i]; }
  G4double GetEmittedGammaPosY(G4int i){ return fEmittedGammaPosY[i]; }
  G4double GetEmittedGammaPosZ(G4int i){ return fEmittedGammaPosZ[i]; }
  G4double GetEmittedGammaPhi(G4int i){ return fEmittedGammaPhi[i]; }
  G4double GetEmittedGammaTheta(G4int i){ return fEmittedGammaTheta[i]; }
  G4int GetNEmittedGammas(){return fNEmittedGammas;}
  G4double GetBeta(G4int i){ return fBeta[i]; }
  G4int GetNBetas(){return fNBetas; }
  G4int GetFullEnergy(){ return fFullEnergy; }
  G4double GetATA(){return fata;}
  G4double GetBTA(){return fbta;}
  G4double GetDTA(){return fdta;}
  G4double GetYTA(){return fyta;}
  G4ThreeVector* GetExitPos(){ return &fExitPos; }
  G4double GetExitBeta(){ return fExitBeta; }
  G4double GetExitTheta(){ return fExitTheta; }
  G4double GetExitPhi(){ return fExitPhi; }
  G4double GetExitTime(){ return fExitTime; }
  G4int    GetFilterCode(){return ffiltercode;}
  G4bool   WriteEvent(){return fwrite;}

private:

  G4double fEmittedGammaEnergies[MAX_SIM_GAMMAS];
  G4double fEmittedGammaPosX[MAX_SIM_GAMMAS];
  G4double fEmittedGammaPosY[MAX_SIM_GAMMAS];
  G4double fEmittedGammaPosZ[MAX_SIM_GAMMAS];
  G4double fEmittedGammaPhi[MAX_SIM_GAMMAS];
  G4double fEmittedGammaTheta[MAX_SIM_GAMMAS];
  G4int    fNEmittedGammas;
  G4double fBeta[MAX_SIM_GAMMAS];
  G4int    fBetaTrackID[MAX_SIM_GAMMAS];
  G4int    fNBetas;
  G4int    fFullEnergy;
  G4double fata;
  G4double fbta;
  G4double fdta;
  G4double fyta;
  G4ThreeVector fExitPos;
  G4double fExitBeta;
  G4double fExitTheta;
  G4double fExitPhi;
  G4double fExitTime;
  G4int    ffiltercode;
  G4bool   fwrite;
};

#endif
