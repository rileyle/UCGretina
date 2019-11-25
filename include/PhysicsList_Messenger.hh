#ifndef PHYSICSLIST_MESSENGER_H
#define PHYSICSLIST_MESSENGER_H

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"

class PhysicsList;

class PhysicsList_Messenger: public G4UImessenger
{
public:
  PhysicsList_Messenger(PhysicsList*);
  ~PhysicsList_Messenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  PhysicsList* aPhysicsList;

  G4UIdirectory* physDir;

  G4UIcmdWithABool* CorrCmd;
  G4UIcmdWithABool* PolCmd;
  G4UIcmdWithAString* AddCmd;
};

#endif//PHYSICSLIST_MESSENGER_H
