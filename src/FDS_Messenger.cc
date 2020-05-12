#ifndef SCANNING
#include "FDS_Messenger.hh"

FDS_Messenger::FDS_Messenger(FDS* pTarget)
:myTarget(pTarget)
{
  const char *aLine;
  G4String commandName;
 
  FDSDir = new G4UIdirectory("/FDS/");
  FDSDir->SetGuidance("FRIB Decay Station control.");
  CloverDir = new G4UIdirectory("/FDS/Clover/");
  CloverDir->SetGuidance("FRIB Decay Station clover array control.");
  ShieldDir = new G4UIdirectory("/FDS/Clover/Shield/");
  ShieldDir->SetGuidance("FRIB Decay Station clover shield control.");
  LaBrDir = new G4UIdirectory("/FDS/LaBr/");
  LaBrDir->SetGuidance("FRIB Decay Station LaBr array control.");
  
  commandName = "/FDS/Clover/EulerFile";
  aLine = commandName.c_str();
  CloverEulerCmd = new G4UIcmdWithAString(aLine, this);
  CloverEulerCmd->SetGuidance("Set the clover Euler-angle file for the FDS.");
  CloverEulerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/FDS/Clover/OuterDeadLayer";
  aLine = commandName.c_str();
  CloverOuterDLCmd = new G4UIcmdWithADoubleAndUnit(aLine, this);
  CloverOuterDLCmd->SetGuidance("Set the clover outer dead layer thickness.");
  CloverOuterDLCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  commandName = "/FDS/Clover/OuterDeadLayer";

  commandName = "/FDS/Clover/CoaxialDeadLayer";
  aLine = commandName.c_str();
  CloverCoaxDLCmd = new G4UIcmdWithADoubleAndUnit(aLine, this);
  CloverCoaxDLCmd->SetGuidance("Set the clover coaxial dead layer thickness.");
  CloverCoaxDLCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/FDS/Clover/Shield/EulerFile";
  aLine = commandName.c_str();
  ShieldEulerCmd = new G4UIcmdWithAString(aLine, this);
  ShieldEulerCmd->SetGuidance("Set the clover shield Euler-angle file for the FDS.");
  ShieldEulerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/FDS/LaBr/EulerFile";
  aLine = commandName.c_str();
  LaBrEulerCmd = new G4UIcmdWithAString(aLine, this);
  LaBrEulerCmd->SetGuidance("Set the LaBr Euler-angle file for the FDS.");
  LaBrEulerCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}



FDS_Messenger::~FDS_Messenger()
{
  delete CloverEulerCmd;
  delete CloverOuterDLCmd;
  delete CloverCoaxDLCmd;
  delete ShieldEulerCmd;
  delete LaBrEulerCmd;
}


void FDS_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == CloverEulerCmd ) {
    myTarget->SetCloverEuler(newValue);
  }

  if( command == CloverOuterDLCmd ) {
    myTarget->SetCloverOuterDL(CloverOuterDLCmd->GetNewDoubleValue(newValue));
  }

  if( command == CloverCoaxDLCmd ) {
    myTarget->SetCloverCoaxialDL(CloverCoaxDLCmd->GetNewDoubleValue(newValue));
  }

  if( command == ShieldEulerCmd ) {
    myTarget->SetShieldEuler(newValue);
  }

  if( command == LaBrEulerCmd ) {
    myTarget->SetLaBrEuler(newValue);
  }

}

#endif
