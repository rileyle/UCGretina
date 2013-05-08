#ifdef LHTARGET
#include "Target_LH_Messenger.hh"

Target_Messenger::Target_Messenger(Target* Tar)
:aTarget(Tar)
{ 
 
  TargetDir = new G4UIdirectory("/Target/");
  TargetDir->SetGuidance("Target control.");
  
  CCmd = new G4UIcmdWithAString("/Target/Cell",this);
  CCmd->SetGuidance("Select the target cell (thick: 200 mg/cm^2, thin: 50 mg/cm^2, empty: cell body, no windows)");
  CCmd->SetParameterName("choice",false);
  CCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  BCmd = new G4UIcmdWithADoubleAndUnit("/Target/Bulge",this);
  BCmd->SetGuidance("Set the thickness of the bulges.)");
  BCmd->SetParameterName("choice",false);
  BCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ACmd = new G4UIcmdWithADoubleAndUnit("/Target/Angle",this);
  ACmd->SetGuidance("Set the angle of tilt of the target.)");
  ACmd->SetParameterName("choice",false);
  ACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MatCmd = new G4UIcmdWithAString("/Target/Material",this);
  MatCmd->SetGuidance("Select material for the target.");
  MatCmd->SetParameterName("choice",false);
  MatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ScDTarCmd = new G4UIcmdWithADouble("/Target/ScaleDensity",this);
  ScDTarCmd->SetGuidance("Scale target density/stopping powers.");
  ScDTarCmd->SetParameterName("choice",false);
  ScDTarCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NSCmd = new G4UIcmdWithAnInteger("/Target/NStep",this);
  NSCmd->SetGuidance("Select the number of steps in the target");
  NSCmd->SetParameterName("choice",false);
  NSCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  sFCmd = new G4UIcmdWithAString("/Target/sourceFrame",this);
  sFCmd->SetGuidance("Select a source frame (eu152_Z2707, cs137_E2879, or co56_2012).");
  sFCmd->SetParameterName("choice",false);
  sFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RepCmd = new G4UIcmdWithoutParameter("/Target/Report",this);
  RepCmd->SetGuidance("Report target parameters");  

  sledCmd = new G4UIcmdWithoutParameter("/Target/Sled",this);
  sledCmd->SetGuidance("Construct the target sled.");  

}



Target_Messenger::~Target_Messenger()
{
  delete TargetDir;
  delete CCmd;
  delete BCmd;
  delete ACmd;
  delete MatCmd;
  delete ScDTarCmd;
  delete NSCmd;
  delete sFCmd;
  delete sledCmd;
  delete RepCmd;
}


void Target_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == CCmd )
    { aTarget->setTargetCell(newValue); }
  if( command == BCmd )
    { aTarget->setBulgeThickness(BCmd->GetNewDoubleValue(newValue)); }
  if( command == ACmd )
    { aTarget->setTargetAngle(ACmd->GetNewDoubleValue(newValue)); }
  if( command == MatCmd )
   { aTarget->setMaterial(newValue);} 
  if( command == ScDTarCmd )
    { aTarget->ScaleDensity(ScDTarCmd->GetNewDoubleValue(newValue)); }
  if( command == NSCmd )
    { aTarget->setNStep(NSCmd->GetNewIntValue(newValue)); }
  if( command == sFCmd )
    { aTarget->setSourceFrame(newValue); }
  if( command == sledCmd )
    { aTarget->setSled(); }
  if( command == RepCmd )
    { aTarget->Report(); }
}
#endif
