#ifndef TrackerIonSD_Messenger_h
#define TrackerIonSD_Messenger_h 1

#include "TrackerIonSD.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"


class TrackerIonSD_Messenger: public G4UImessenger
{

  public:
    TrackerIonSD_Messenger(TrackerIonSD*);
   ~TrackerIonSD_Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    TrackerIonSD* tracker;    
    G4UIdirectory*             PrtIDir;
    G4UIcmdWithoutParameter*   PrtISCmd;
    G4UIcmdWithoutParameter*   PrtIUCmd;

};


#endif

