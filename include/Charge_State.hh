#ifndef Charge_State_h
#define Charge_State_h 1

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

class Charge_State
{
public:
  Charge_State();
  ~Charge_State();
  void  SetCharge(G4int q){Charge=q;}
  void  SetUnReactedFraction(G4double f){UnReactedFraction=f;}
  void  SetReactedFraction(G4double f){ReactedFraction=f;}
  void  SetKE(G4double e){setKE=e;useSetKEu=false;}
  void  SetKEu(G4double e,G4int A){setKEu=e;useSetKEu=true;setKE=setKEu*A;}
  G4int GetCharge(){return Charge;}
  G4double GetUnReactedFraction(){return UnReactedFraction;}
  G4double GetReactedFraction(){return ReactedFraction;}
  G4double GetSetKEu(){return setKEu;}
  G4double GetSetKE(){return setKE;}
  G4bool   GetUseSetKEu(){return useSetKEu;}
private:
  G4int    Charge;
  G4double UnReactedFraction;
  G4double ReactedFraction;
  G4double setKEu;
  G4double setKE;
  G4bool   useSetKEu;
};
#endif


           
