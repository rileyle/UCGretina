#ifndef Incoming_Beam_h
#define Incoming_Beam_h 1

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "Randomize.hh"
class Incoming_Beam
{
public:
  Incoming_Beam();
  ~Incoming_Beam();
  
public:

  void Report();
  void setA(G4int);
  void setZ(G4int);
  void setKE(G4double);
  void setKEu(G4double); 
  void setDTAFile(G4String);
  void setfcX(G4double);
  void setfcDX(G4double);
  void setfcY(G4double);
  void setfcDY(G4double);
  void setfcZ(G4double);
  void setDpp(G4double);
  void setmaxAta(G4double);
  void setmaxBta(G4double);
  void setAta0(G4double a){ata0=a;}
  void setBta0(G4double b){bta0=b;}


  G4int getA()    {return A;}
  G4int getZ()    {return Z;}
  G4double getEx(){return Ex;}
  G4double getKE(G4ParticleDefinition*);
  G4double getAta0(){return ata0;}
  G4double getBta0(){return bta0;}
  G4ThreeVector getPosition();
  G4ThreeVector getDirection();
private:
  
  G4int A;
  G4int Z;
  G4double Ex;
  G4double KE;
  G4double KEu;
  G4String dtaFileName;
  G4double dta[1000];
  G4int Ndta;
  G4double dtaMin;
  G4double dtaMax;
  G4double dtaBin;
  G4double Dpp;
  G4double fcX;
  G4double fcDX;
  G4double fcY;
  G4double fcDY;
  G4double fcZ;
  G4double maxAta;
  G4double maxBta;
  G4double ata0;
  G4double bta0;

};

#endif
