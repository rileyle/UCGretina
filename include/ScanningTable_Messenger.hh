#ifdef SCANNING
#ifndef Scanning_Table_Messenger_h
#define Scanning_Table_Messenger_h

#include "ScanningTable.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

class ScanningTable_Messenger: public G4UImessenger
{
  public:
  ScanningTable_Messenger(ScanningTable*);
  ~ScanningTable_Messenger();
    
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  ScanningTable* scanningTable;

  G4UIdirectory*  ScanningTableDir;
  G4UIcmdWithAString* CADPathCmd;
  G4UIcmdWithoutParameter* ConstructCmd;
  G4UIcmdWithoutParameter* CloverCartCmd;
  G4UIcmdWithoutParameter* CartFrameCmd;  
  G4UIcmdWithoutParameter* SlitCmd;
  G4UIcmdWithoutParameter* SlitMountCmd;  
  G4UIcmdWithoutParameter* CollimatorCmd;
  G4UIcmdWithoutParameter* CollimatorInsertCmd;
  G4UIcmdWithoutParameter* CollimatorMountCmd;
  G4UIcmdWithoutParameter* CuTargetCmd;
  G4UIcmdWithoutParameter* ShieldCmd;
  G4UIcmdWithADoubleAndUnit* ControllerXCmd;
  G4UIcmdWithADoubleAndUnit* ControllerYCmd;
  G4UIcmdWithADoubleAndUnit* ControllerZCmd;
  G4UIcmdWithADoubleAndUnit* CloverZCmd;
  G4UIcmdWithADoubleAndUnit* CollRCmd;
  G4UIcmdWithADoubleAndUnit* SlitWCmd;
};

#endif
#endif

