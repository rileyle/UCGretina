#include "TrackerIonSD_Messenger.hh"


TrackerIonSD_Messenger::TrackerIonSD_Messenger(TrackerIonSD* TISD)
:tracker(TISD) 
{ 
 
  PrtIDir = new G4UIdirectory("/IonPrint/");
  PrtIDir->SetGuidance("Event by event print control for.");
  
  PrtISCmd = new G4UIcmdWithoutParameter("/IonPrint/Track_Set",this);
  PrtISCmd->SetGuidance("Sets printing of track ion results.");

  PrtIUCmd = new G4UIcmdWithoutParameter("/IonPrint/Track_UnSet",this);
  PrtIUCmd->SetGuidance("Un sets printing of track ion results.");
  
}



TrackerIonSD_Messenger::~TrackerIonSD_Messenger()
{
  
  delete PrtIDir;
  delete PrtISCmd;
  delete PrtIUCmd;

}



void TrackerIonSD_Messenger::SetNewValue(G4UIcommand* command,G4String)
{ 


  if( command == PrtISCmd )
    {tracker ->SetPrint();}

  if( command == PrtIUCmd )
    {tracker ->UnSetPrint();}

}

