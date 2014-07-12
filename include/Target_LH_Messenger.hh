#ifndef Target_Messenger_h
#define Target_Messenger_h 1

#include "Target_LH.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"

class Target_Messenger: public G4UImessenger
{
public:
  Target_Messenger(Target*);
  ~Target_Messenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  Target* aTarget;

  G4UIdirectory*             TargetDir;
  G4UIcmdWithAString*        CCmd; 
  G4UIcmdWithADoubleAndUnit* BCmd;
  G4UIcmdWithoutParameter*   WCmd; 
  G4UIcmdWithADoubleAndUnit* ACmd;
  G4UIcmdWithAString*        MatCmd;  
  G4UIcmdWithADouble*        SDTarCmd;
  G4UIcmdWithAnInteger*      NSCmd;
  G4UIcmdWithAString*        sFCmd; 
  G4UIcmdWithoutParameter*   sledCmd; 
  G4UIcmdWithoutParameter*   RepCmd; 
  G4UIcmdWithoutParameter*   GCmd;
};

#endif

