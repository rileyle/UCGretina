#ifndef Incoming_Beam_Messenger_h
#define Incoming_Beam_Messenger_h 1

#include "Incoming_Beam.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"


class Incoming_Beam_Messenger: public G4UImessenger
{
  public:
    Incoming_Beam_Messenger(Incoming_Beam*);
   ~Incoming_Beam_Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Incoming_Beam* BeamIn;    
    G4UIdirectory*             BeamInDir;
    G4UIdirectory*             fcDir;
    G4UIcmdWithAnInteger*      ACmd;
    G4UIcmdWithAnInteger*      ZCmd;
    G4UIcmdWithADoubleAndUnit* ExCmd;
    G4UIcmdWithADoubleAndUnit* KECmd;
    G4UIcmdWithADoubleAndUnit* KEuCmd;
    G4UIcmdWithAString*        dtaCmd;
    G4UIcmdWithAString*        pdCmd;
    G4UIcmdWithADoubleAndUnit* fcXCmd;
    G4UIcmdWithADoubleAndUnit* fcDXCmd;
    G4UIcmdWithADoubleAndUnit* fcYCmd;
    G4UIcmdWithADoubleAndUnit* fcDYCmd;
    G4UIcmdWithADoubleAndUnit* fcZCmd;
    G4UIcmdWithADoubleAndUnit* maxACmd;
    G4UIcmdWithADoubleAndUnit* maxBCmd;
    G4UIcmdWithADoubleAndUnit* Ata0Cmd;
    G4UIcmdWithADoubleAndUnit* Bta0Cmd;
    G4UIcmdWithADouble*        DppCmd;
    G4UIcmdWithoutParameter*   RepCmd;
  

};


#endif

