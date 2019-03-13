#include "Reaction_Messenger.hh"
#include "Reaction.hh"

Reaction_Messenger::Reaction_Messenger(Reaction* rxn) : theReaction(rxn){
  reactionDir = new G4UIdirectory("/reaction/");
  reactionDir->SetGuidance("Reaction control");

  popCmd = new G4UIcommand("/reaction/population",this);
  popCmd->SetGuidance("Set the population fraction of a magnetic substate");
  popCmd->SetGuidance("[usage] /reaction/population 2M(int) frac(double)");
  G4UIparameter* param;
  param = new G4UIparameter("twoM",'i',false);
  popCmd->SetParameter(param);
  param = new G4UIparameter("frac",'d',false);
  popCmd->SetParameter(param);
}

Reaction_Messenger::~Reaction_Messenger(){
  delete reactionDir;
  delete popCmd;
}

void Reaction_Messenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
  if(cmd == popCmd){
    G4Tokenizer next(newValue);
    G4int ftwoM = StoI(next());
    G4double ffrac = StoD(next());
    theReaction->SetPopulation(ftwoM,ffrac);
  }
}
