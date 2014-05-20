#ifndef LHTARGET
#include "Beam_Tube_Messenger.hh"

Beam_Tube_Messenger::Beam_Tube_Messenger(Beam_Tube* BT)
:BeamTube(BT)
{ 
 
  BeamTubeDir = new G4UIdirectory("/BeamTube/");
  BeamTubeDir->SetGuidance("Beam tube control.");
  
  MatCmd = new G4UIcmdWithAString("/BeamTube/Material",this);
  MatCmd->SetGuidance("Select material for the beam tube.");
  MatCmd->SetParameterName("choice",false);
  MatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  RminCmd = new G4UIcmdWithADoubleAndUnit("/BeamTube/R_min",this);
  RminCmd->SetGuidance("Select the inner radius of the beam tube");
  RminCmd->SetParameterName("choice",false);
  RminCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RmaxCmd = new G4UIcmdWithADoubleAndUnit("/BeamTube/R_max",this);
  RmaxCmd->SetGuidance("Select the outer radius of the beam tube");
  RmaxCmd->SetParameterName("choice",false);
  RmaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  LengthCmd = new G4UIcmdWithADoubleAndUnit("/BeamTube/Length",this);
  LengthCmd->SetGuidance("Select the length of the beam tube");
  LengthCmd->SetParameterName("choice",false);
  LengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RepCmd = new G4UIcmdWithoutParameter("/BeamTube/Report",this);
  RepCmd->SetGuidance("Report beam tube parameters");  

}



Beam_Tube_Messenger::~Beam_Tube_Messenger()
{
  delete LengthCmd;
  delete RminCmd;
  delete RmaxCmd;
  delete MatCmd;
  delete BeamTubeDir;
  delete RepCmd;
}


void Beam_Tube_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == MatCmd )
   { BeamTube->setMaterial(newValue);} 
  if( command == RminCmd )
   { BeamTube->setRmin(RminCmd->GetNewDoubleValue(newValue));}
  if( command == RmaxCmd )
   { BeamTube->setRmax(RmaxCmd->GetNewDoubleValue(newValue));}
  if( command == LengthCmd )
   { BeamTube->setLength(LengthCmd->GetNewDoubleValue(newValue));}
  if( command == RepCmd )
   { BeamTube->Report();}

}
#endif
