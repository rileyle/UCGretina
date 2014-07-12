#ifndef LHTARGET
#ifndef Greta_Chamber_Messenger_h
#define Greta_Chamber_Messenger_h 1

#include "Greta_Chamber.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

class Greta_Chamber_Messenger: public G4UImessenger
{
  public:
    Greta_Chamber_Messenger(Greta_Chamber*);
   ~Greta_Chamber_Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Greta_Chamber* GretaChamber;
    
    G4UIdirectory*             GretaChamberDir;
    G4UIcmdWithADoubleAndUnit* RmaxCmd;
    G4UIcmdWithADoubleAndUnit* RminCmd;
    G4UIcmdWithoutParameter*   RepCmd;

};

#endif
#endif

