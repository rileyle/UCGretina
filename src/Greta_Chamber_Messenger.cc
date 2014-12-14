#ifndef LHTARGET
#include "Greta_Chamber_Messenger.hh"

Greta_Chamber_Messenger::Greta_Chamber_Messenger(Greta_Chamber* GC)
:GretaChamber(GC)
{ 
 
  GretaChamberDir = new G4UIdirectory("/GretaChamber/");
  GretaChamberDir->SetGuidance("Greta chamber control.");

  RminCmd = new G4UIcmdWithADoubleAndUnit("/GretaChamber/R_min",this);
  RminCmd->SetGuidance("Select the inner radius of the GRETA chamber");
  RminCmd->SetParameterName("choice",false);
  RminCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RmaxCmd = new G4UIcmdWithADoubleAndUnit("/GretaChamber/R_max",this);
  RmaxCmd->SetGuidance("Select the outer radius of the GRETA chamber");
  RmaxCmd->SetParameterName("choice",false);
  RmaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  CutCmd = new G4UIcmdWithoutParameter("/GretaChamber/Cutaway",this);
  CutCmd->SetGuidance("Build a cutaway view of the GRETA chamber. For visualization only!");  

  RepCmd = new G4UIcmdWithoutParameter("/GretaChamber/Report",this);
  RepCmd->SetGuidance("Report beam tube parameters");  

}



Greta_Chamber_Messenger::~Greta_Chamber_Messenger()
{
  delete GretaChamberDir;
  delete RminCmd;
  delete RmaxCmd;
  delete CutCmd;
  delete RepCmd;
}


void Greta_Chamber_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == RminCmd )
   { GretaChamber->setRmin(RminCmd->GetNewDoubleValue(newValue));}
  if( command == RmaxCmd )
   { GretaChamber->setRmax(RmaxCmd->GetNewDoubleValue(newValue));}
  if( command == CutCmd )
   { GretaChamber->setCutaway();}
  if( command == RepCmd )
   { GretaChamber->Report();}
}
#endif
