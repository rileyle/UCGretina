#ifndef SCANNING
#ifndef FDS_Messenger_h
#define FDS_Messenger_h 1

#include "FDS.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

class FDS;

class FDS_Messenger: public G4UImessenger
{
public:
  FDS_Messenger(FDS*);
  ~FDS_Messenger();
    
private:
  FDS* myTarget;

  G4UIdirectory* FDSDir;
  G4UIdirectory* CloverDir;
  G4UIdirectory* ShieldDir;
  G4UIdirectory* LaBrDir;

  G4UIcmdWithAString* CloverEulerCmd;
  G4UIcmdWithADoubleAndUnit* CloverOuterDLCmd;
  G4UIcmdWithADoubleAndUnit* CloverCoaxDLCmd;
  G4UIcmdWithAString* ShieldEulerCmd;
  G4UIcmdWithAString* LaBrEulerCmd;

public:
  void SetNewValue(G4UIcommand*, G4String);
};

#endif
#endif
