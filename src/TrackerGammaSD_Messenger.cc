#include "TrackerGammaSD_Messenger.hh"


TrackerGammaSD_Messenger::TrackerGammaSD_Messenger(TrackerGammaSD* TGSD)
:tracker(TGSD) 
{ 

  PrtGDir = new G4UIdirectory("/GammaPrint/");
  PrtGDir->SetGuidance("Event by event print control for.");

  PrtGSCmd = new G4UIcmdWithoutParameter("/GammaPrint/Track_Set",this);
  PrtGSCmd->SetGuidance("Sets printing of track gamma results.");

  PrtGUCmd = new G4UIcmdWithoutParameter("/GammaPrint/Track_UnSet",this);
  PrtGUCmd->SetGuidance("Un sets printing of track gamma results.");

  MethodCmd = new G4UIcmdWithAnInteger("/Tracking/Method",this);
  MethodCmd->SetGuidance("Set the tracking method: 0 to track the primary and secondary partices, 1 to track only the primary gammas (and gammas from annihilation events).");
}



TrackerGammaSD_Messenger::~TrackerGammaSD_Messenger()

{  

  delete PrtGDir;
  delete PrtGSCmd;
  delete PrtGUCmd;
  delete MethodCmd;

}



void TrackerGammaSD_Messenger::SetNewValue(G4UIcommand* command,
					   G4String newValue)
{ 

  if( command == PrtGSCmd )
    { tracker->SetPrint(); }

  if( command == PrtGUCmd )
    { tracker->UnSetPrint(); }

  if( command == MethodCmd )
    { tracker->SetTrackingMethod(MethodCmd->GetNewIntValue(newValue)); }

}

