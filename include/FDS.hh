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
    
  void SetFDSCloverEuler(G4String file){fdsCloverEulerFile = file;}
  void SetFDSShieldEuler(G4String file){fdsShieldEulerFile = file;}
  void SetFDSLaBrEuler(G4String file){fdsLaBrEulerFile = file;}

private:
  FDS_Messenger *myMessenger;

  void ReadFDSCloverEulerFile();
  void ReadFDSShieldEulerFile();
  void ReadFDSLaBrEulerFile();

  Materials* materials;
  
  G4int nCloverEuler;
  G4int nShieldEuler;
  G4int nLaBrEuler;
  std::vector<CeulerAngles> cloverEuler;
  std::vector<CeulerAngles> shieldEuler;
  std::vector<CeulerAngles> labrEuler;
  
  G4String fdsCloverEulerFile;
  G4String fdsShieldEulerFile;
  G4String fdsLaBrEulerFile;

};

#endif
#endif
