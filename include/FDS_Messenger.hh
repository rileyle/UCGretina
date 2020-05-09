#ifndef SCANNING
#ifndef FDS_Messenger_h
#define FDS_Messenger_h 1

#include "FDS.hh"

#include "G4UIcmdWithAString.hh"

class FDS;

class FDS_Messenger: public G4UImessenger
{
public:
  FDS_Messenger(FDS*);
  ~FDS_Messenger();
    
private:
  FDS* myTarget;

  G4UIdirectory* FDSDir;
  G4UIdirectory* FDSCloverDir;
  G4UIdirectory* FDSShieldDir;
  G4UIdirectory* FDSLaBrDir;

  G4UIcmdWithAString* FDSCloverEulerCmd;
  G4UIcmdWithAString* FDSShieldEulerCmd;
  G4UIcmdWithAString* FDSLaBrEulerCmd;

public:
  void SetNewValue(G4UIcommand*, G4String);
};

#endif
#endif
