#ifndef LHTARGET
#include "Target_Messenger.hh"

Target_Messenger::Target_Messenger(Target* Tar)
:aTarget(Tar)
{ 
 
  TargetDir = new G4UIdirectory("/Target/");
  TargetDir->SetGuidance("Target control.");
  
  MatCmd = new G4UIcmdWithAString("/Target/Material",this);
  MatCmd->SetGuidance("Select material for the target.");
  MatCmd->SetParameterName("choice",false);
  MatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  XCmd = new G4UIcmdWithADoubleAndUnit("/Target/X_length",this);
  XCmd->SetGuidance("Select the X length for the target");
  XCmd->SetParameterName("choice",false);
  XCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  YCmd = new G4UIcmdWithADoubleAndUnit("/Target/Y_length",this);
  YCmd->SetGuidance("Select the Y length for the target");
  YCmd->SetParameterName("choice",false);
  YCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ZCmd = new G4UIcmdWithADoubleAndUnit("/Target/Thickness",this);
  ZCmd->SetGuidance("Select the thickness for the target");
  ZCmd->SetParameterName("choice",false);
  ZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ScDTarCmd = new G4UIcmdWithADouble("/Target/ScaleDensity",this);
  ScDTarCmd->SetGuidance("Scale target density/stopping powers.");
  ScDTarCmd->SetParameterName("choice",false);
  ScDTarCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PosXCmd = new G4UIcmdWithADoubleAndUnit("/Target/SetPosition_X",this);
  PosXCmd->SetGuidance("Select the position of the target along the X axis.");
  PosXCmd->SetParameterName("choice",false);
  PosXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PosYCmd = new G4UIcmdWithADoubleAndUnit("/Target/SetPosition_Y",this);
  PosYCmd->SetGuidance("Select the position of the target along the Y axis.");
  PosYCmd->SetParameterName("choice",false);
  PosYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PosZCmd = new G4UIcmdWithADoubleAndUnit("/Target/SetPosition_Z",this);
  PosZCmd->SetGuidance("Select the position of the target along the beam axis (Z direction.");
  PosZCmd->SetParameterName("choice",false);
  PosZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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
  delete XCmd;
  delete YCmd;
  delete ZCmd;
  delete MatCmd;
  delete TargetDir;
  delete RepCmd;
  delete PosXCmd;
  delete PosYCmd;
  delete PosZCmd;
  delete ScDTarCmd;
  delete NSCmd;
  delete sFCmd;
  delete sledCmd;
}


void Target_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == MatCmd )
   { aTarget->setMaterial(newValue);} 
  if( command == XCmd )
   { aTarget->setX(XCmd->GetNewDoubleValue(newValue));}
  if( command == YCmd )
   { aTarget->setY(YCmd->GetNewDoubleValue(newValue));}
  if( command == ZCmd )
   { aTarget->setZ(ZCmd->GetNewDoubleValue(newValue));}
  if( command == ScDTarCmd )
   { aTarget->ScaleDensity(ScDTarCmd->GetNewDoubleValue(newValue));}
  if( command == PosXCmd )
   { aTarget->SetPosition(PosXCmd->GetNewDoubleValue(newValue), 0., 0.);}
  if( command == PosYCmd )
   { aTarget->SetPosition(0., PosYCmd->GetNewDoubleValue(newValue), 0.);}
  if( command == PosZCmd )
   { aTarget->SetPosition(0., 0., PosZCmd->GetNewDoubleValue(newValue));}
  if( command == NSCmd )
   { aTarget->setNStep(NSCmd->GetNewIntValue(newValue));}
  if( command == RepCmd )
   { aTarget->Report();}
  if( command == sFCmd )
   { aTarget->setSourceFrame(newValue);} 
  if( command == sledCmd )
   { aTarget->setSled();} 

}
#endif
