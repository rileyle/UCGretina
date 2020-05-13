#ifndef SCANNING
#ifndef FDS_h
#define FDS_h 1

#include "G4VUserDetectorConstruction.hh"
#include "DetectorConstruction.hh"
#include "TrackerGammaSD.hh"
#include "RunAction.hh"

#include "G4SystemOfUnits.hh"

#include "Gretina_Helper.hh"
#include "Clover_Detector.hh"
#include "FDS_Messenger.hh"

using namespace std;

class FDS_Messenger;

class FDS
{
public:
  FDS();
  ~FDS();

  void Placement(G4LogicalVolume*);
    
  void SetCloverEuler(G4String file){cloverEulerFile = file;}
  void SetCloverOuterDL(G4double t){cloverOuterDL = t;}
  void SetCloverCoaxialDL(G4double t){cloverCoaxialDL = t;}
  void SetShieldEuler(G4String file){shieldEulerFile = file;}
  void SetLaBrEuler(G4String file){labrEulerFile = file;}

private:
  FDS_Messenger *myMessenger;

  void ReadCloverEulerFile();
  void ReadShieldEulerFile();
  void ReadLaBrEulerFile();

  G4int nCloverEuler;
  G4int nShieldEuler;
  G4int nLaBrEuler;
  std::vector<CeulerAngles> cloverEuler;
  std::vector<CeulerAngles> shieldEuler;
  std::vector<CeulerAngles> labrEuler;
  
  G4String cloverEulerFile;
  G4String shieldEulerFile;
  G4String labrEulerFile;

  G4double cloverOuterDL;
  G4double cloverCoaxialDL;

};

#endif
#endif
