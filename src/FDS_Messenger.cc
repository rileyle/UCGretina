#ifndef SCANNING
#include "FDS_Messenger.hh"

FDS_Messenger::FDS_Messenger(FDS* pTarget)
:myTarget(pTarget)
{
  const char *aLine;
  G4String commandName;
 
  FDSDir = new G4UIdirectory("/FDS/");
  FDSDir->SetGuidance("FRIB Decay Station control.");
  FDSCloverDir = new G4UIdirectory("/FDS/Clover/");
  FDSCloverDir->SetGuidance("FRIB Decay Station clover array control.");
  FDSShieldDir = new G4UIdirectory("/FDS/Clover/Shield/");
  FDSShieldDir->SetGuidance("FRIB Decay Station clover shield control.");
  FDSLaBrDir = new G4UIdirectory("/FDS/LaBr/");
  FDSLaBrDir->SetGuidance("FRIB Decay Station LaBr array control.");
  
  commandName = "/FDS/Clover/EulerFile";
  aLine = commandName.c_str();
  FDSCloverEulerCmd = new G4UIcmdWithAString(aLine, this);
  FDSCloverEulerCmd->SetGuidance("Set the clover Euler-angle file for the FDS.");
  FDSCloverEulerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/FDS/Clover/Shield/EulerFile";
  aLine = commandName.c_str();
  FDSShieldEulerCmd = new G4UIcmdWithAString(aLine, this);
  FDSShieldEulerCmd->SetGuidance("Set the clover shield Euler-angle file for the FDS.");
  FDSShieldEulerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/FDS/LaBr/EulerFile";
  aLine = commandName.c_str();
  FDSLaBrEulerCmd = new G4UIcmdWithAString(aLine, this);
  FDSLaBrEulerCmd->SetGuidance("Set the LaBr Euler-angle file for the FDS.");
  FDSLaBrEulerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}



FDS_Messenger::~FDS_Messenger()
{
  delete FDSCloverEulerCmd;
  delete FDSShieldEulerCmd;
  delete FDSLaBrEulerCmd;
}


void FDS_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == FDSCloverEulerCmd ) {
    myTarget->SetFDSCloverEuler(newValue);
  }

  if( command == FDSShieldEulerCmd ) {
    myTarget->SetFDSShieldEuler(newValue);
  }

  if( command == FDSLaBrEulerCmd ) {
    myTarget->SetFDSLaBrEuler(newValue);
  }

}

#endif
