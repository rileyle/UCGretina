#ifndef REACTION_MESSENGER_H
#define REACTION_MESSENGER_H

#include "G4UImessenger.hh"
#include "G4Tokenizer.hh"

class Reaction;

class Reaction_Messenger : public G4UImessenger{
public:
  Reaction_Messenger(Reaction*);
  ~Reaction_Messenger();

  void SetNewValue(G4UIcommand*,G4String);

private:
  Reaction* theReaction;
  G4UIdirectory* reactionDir;
  G4UIcommand* popCmd;
};

#endif//REACTION_MESSENGER_H
