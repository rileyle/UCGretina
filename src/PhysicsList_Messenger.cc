#include "PhysicsList_Messenger.hh"

#include "PhysicsList.hh"

PhysicsList_Messenger::PhysicsList_Messenger(PhysicsList* pl):aPhysicsList(pl){

  physDir = new G4UIdirectory("/PhysicsList/");
  physDir->SetGuidance("PhysicsList control");

  CorrCmd = new G4UIcmdWithABool("/PhysicsList/AngularCorrelations",this);
  CorrCmd->SetGuidance("Turn gamma-ray angular correlations on/off");
  CorrCmd->SetParameterName("choice",false);
  CorrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}

PhysicsList_Messenger::~PhysicsList_Messenger(){
  delete CorrCmd;
  delete physDir;
}

void PhysicsList_Messenger::SetNewValue(G4UIcommand* command, G4String newValue){
  if( command == CorrCmd){
    aPhysicsList->SetGammaAngularCorrelations(CorrCmd->GetNewBoolValue(newValue));
  }
}
