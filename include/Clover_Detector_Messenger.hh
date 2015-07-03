#ifndef Clover_Detector_Messenger_h
#define Clover_Detector_Messenger_h 1

#include "Clover_Detector.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"

class Clover_Detector_Messenger: public G4UImessenger
{
  public:
    Clover_Detector_Messenger(Clover_Detector*);
   ~Clover_Detector_Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Clover_Detector* CloverDet;
   
    G4UIdirectory*             CloverDir;  

    G4UIcmdWithADoubleAndUnit* XCmd;
    G4UIcmdWithADoubleAndUnit* YCmd;
    G4UIcmdWithADoubleAndUnit* ZCmd;
};


#endif

