#include "Background_Sphere_Messenger.hh"

Background_Sphere_Messenger::Background_Sphere_Messenger(Background_Sphere* BS)
:BackgroundSphere(BS)
{ 
 
  BackgroundSphereDir = new G4UIdirectory("/BackgroundSphere/");
  BackgroundSphereDir->SetGuidance("Background sphere control.");
  
  MatCmd = new G4UIcmdWithAString("/BackgroundSphere/Material",this);
  MatCmd->SetGuidance("Set material for the background sphere.");
  MatCmd->SetParameterName("choice",false);
  MatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  RminCmd = new G4UIcmdWithADoubleAndUnit("/BackgroundSphere/R_min",this);
  RminCmd->SetGuidance("Set the inner radius of the beam tube");
  RminCmd->SetParameterName("choice",false);
  RminCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RmaxCmd = new G4UIcmdWithADoubleAndUnit("/BackgroundSphere/R_max",this);
  RmaxCmd->SetGuidance("Set the outer radius of the beam tube");
  RmaxCmd->SetParameterName("choice",false);
  RmaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RepCmd = new G4UIcmdWithoutParameter("/BackgroundSphere/Report",this);
  RepCmd->SetGuidance("Report background sphere parameters");  

}



Background_Sphere_Messenger::~Background_Sphere_Messenger()
{
  delete RminCmd;
  delete RmaxCmd;
  delete MatCmd;
  delete BackgroundSphereDir;
  delete RepCmd;
}


void Background_Sphere_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == MatCmd )
   { BackgroundSphere->setMaterial(newValue);} 
  if( command == RminCmd )
   { BackgroundSphere->setRmin(RminCmd->GetNewDoubleValue(newValue));}
  if( command == RmaxCmd )
   { BackgroundSphere->setRmax(RmaxCmd->GetNewDoubleValue(newValue));}
  if( command == RepCmd )
   { BackgroundSphere->Report();}

}

