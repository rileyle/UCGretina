#ifndef EventAction_Messenger_h
#define EventAction_Messenger_h 1

#include "EventAction.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"


class EventAction_Messenger: public G4UImessenger
{

  public:
    EventAction_Messenger(EventAction*);
   ~EventAction_Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    EventAction*               theEventAction;  
    G4UIdirectory*             OutputDir;
    G4UIcmdWithAString*        OutFileCmd;
    G4UIcmdWithoutParameter*   OutDetCmd;
    G4UIcmdWithADoubleAndUnit* PosResCmd;
    G4UIcmdWithADoubleAndUnit* ThreshECmd;
    G4UIcmdWithADoubleAndUnit* ThreshDECmd;
    G4UIdirectory*             Mode2Dir;
    G4UIcmdWithAString*        Mode2FileCmd;
    G4UIcmdWithoutParameter*   TimeSortCmd;
    G4UIcmdWithAString*        crmatCmd;
    G4UIcmdWithoutParameter*   crysCmd;
    G4UIcmdWithoutParameter*   coordsCmd;
    G4UIcmdWithADoubleAndUnit* PackResCmd;
    G4UIcmdWithADoubleAndUnit* S800KECmd;
    G4UIcmdWithoutParameter*   allS800Cmd;
    G4UIcmdWithoutParameter*   printCmd;
};


#endif

