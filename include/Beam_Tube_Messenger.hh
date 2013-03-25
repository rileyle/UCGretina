#ifndef LHTARGET
#ifndef Beam_Tube_Messenger_h
#define Beam_Tube_Messenger_h 1

#include "Beam_Tube.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

class Beam_Tube_Messenger: public G4UImessenger
{
  public:
    Beam_Tube_Messenger(Beam_Tube*);
   ~Beam_Tube_Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Beam_Tube* BeamTube;
    
    G4UIdirectory*             BeamTubeDir;
    G4UIcmdWithAString*        MatCmd;  
    G4UIcmdWithADoubleAndUnit* RminCmd;
    G4UIcmdWithADoubleAndUnit* RmaxCmd;
    G4UIcmdWithADoubleAndUnit* LengthCmd;
    G4UIcmdWithoutParameter*   RepCmd;

};

#endif
#endif

