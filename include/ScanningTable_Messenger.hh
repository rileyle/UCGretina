#ifndef LHTARGET
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
  G4UIcmdWithoutParameter* ShieldCmd;
  G4UIcmdWithADoubleAndUnit* XShiftCmd;
  G4UIcmdWithADoubleAndUnit* YShiftCmd;
  G4UIcmdWithADoubleAndUnit* ZShiftCmd;
};

#endif
#endif

