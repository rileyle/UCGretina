#ifndef Target_Messenger_h
#define Target_Messenger_h 1

#ifdef LHTARGET
  #include "Target_LH.hh"
#else
  #include "Target.hh"
#endif
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
  G4UIcmdWithAString*        MatCmd;  
  G4UIcmdWithADoubleAndUnit* XCmd;
  G4UIcmdWithADoubleAndUnit* YCmd;
  G4UIcmdWithADoubleAndUnit* ZCmd;
  G4UIcmdWithADoubleAndUnit* PosZCmd;
  G4UIcmdWithADouble*        ScDTarCmd;
  G4UIcmdWithoutParameter*   RepCmd;
  G4UIcmdWithAnInteger*      NSCmd;
  G4UIcmdWithAString*        sFCmd; 
  G4UIcmdWithoutParameter*   sledCmd; 

};


#endif

