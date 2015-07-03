#ifdef SCANNING
#include "Clover_Detector_Messenger.hh"

Clover_Detector_Messenger::Clover_Detector_Messenger(Clover_Detector* SD)
:CloverDet(SD)
{ 
 
  CloverDir = new G4UIdirectory("/Clover/");
  CloverDir->SetGuidance("Clover control.");

  XCmd = new G4UIcmdWithADoubleAndUnit("/Clover/setX",this);
  XCmd->SetGuidance("Set the x position of the detector (crystal center)");
  XCmd->SetParameterName("choice",false);
  XCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  YCmd = new G4UIcmdWithADoubleAndUnit("/Clover/setY",this);
  YCmd->SetGuidance("Set the y position of the detector (crystal center)");
  YCmd->SetParameterName("choice",false);
  YCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ZCmd = new G4UIcmdWithADoubleAndUnit("/Clover/setZ",this);
  ZCmd->SetGuidance("Set the z position of the detector (crystal center)");
  ZCmd->SetParameterName("choice",false);
  ZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}



Clover_Detector_Messenger::~Clover_Detector_Messenger()
{
  delete CloverDir;
  delete XCmd;
  delete YCmd;
  delete ZCmd;
}


void Clover_Detector_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == XCmd )
    {CloverDet->setX(XCmd->GetNewDoubleValue(newValue));}
  if( command == YCmd )
    {CloverDet->setY(YCmd->GetNewDoubleValue(newValue));}
  if( command == ZCmd )
    {CloverDet->setZ(ZCmd->GetNewDoubleValue(newValue));}
}

#endif
