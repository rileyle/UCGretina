#include "EventAction_Messenger.hh"


EventAction_Messenger::EventAction_Messenger(EventAction* EA):theEventAction(EA) 
{ 

  OutputDir = new G4UIdirectory("/Output/");
  OutputDir->SetGuidance("Event by event output control.");

  OutFileCmd = new G4UIcmdWithAString("/Output/Filename",this);
  OutFileCmd->SetGuidance("Output file name.");

  Mode2Dir = new G4UIdirectory("/Mode2/");
  Mode2Dir->SetGuidance("Parameters for simulating Mode2 data.");

  Mode2FileCmd = new G4UIcmdWithAString("/Mode2/Filename",this);
  Mode2FileCmd->SetGuidance("Mode 2 output file name");

  crmatCmd = new G4UIcmdWithAString("/Mode2/crmatFile",this);
  crmatCmd->SetGuidance("Transformations from crystal to world coordinates for Mode2 data (expected in crystal coordinates). Leave unset to get interaction points in world coordinates.");

  coordsCmd = new G4UIcmdWithoutParameter("/Mode2/GretinaCoords",this);
  coordsCmd->SetGuidance("Write interaction points in GRETINA coordinate system --- x down, z beam (standard for Mode2 data).");

  PosResCmd = new G4UIcmdWithADoubleAndUnit("/Mode2/PositionRes",this);
  PosResCmd->SetGuidance("Position resolution for simulating Mode2 data.");

  S800KECmd = new G4UIcmdWithADoubleAndUnit("/Mode2/S800KE",this);
  S800KECmd->SetGuidance("Average KE for computing dta for writing S800 events.");

  printCmd = new G4UIcmdWithoutParameter("/Mode2/Print",this);
  printCmd->SetGuidance("Write mode 2 event information to stdout.");

}



EventAction_Messenger::~EventAction_Messenger()

{  

  delete OutputDir;
  delete OutFileCmd;
  delete Mode2Dir;
  delete Mode2FileCmd;
  delete crmatCmd;
  delete coordsCmd;
  delete PosResCmd;
  delete S800KECmd;
  delete printCmd;

}



void EventAction_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == OutFileCmd )
    {theEventAction->SetOutFile(newValue);}
  if( command == Mode2FileCmd )
    {theEventAction->SetMode2File(newValue);}
  if( command == crmatCmd )
    {theEventAction->SetCrmatFile(newValue);}
  if( command == coordsCmd )
    {theEventAction->SetGretinaCoords();}
  if( command == PosResCmd )
    {theEventAction->SetPosRes(PosResCmd->GetNewDoubleValue(newValue));}
  if( command == S800KECmd )
    {theEventAction->SetS800KE(PosResCmd->GetNewDoubleValue(newValue));}
  if( command == printCmd )
    {theEventAction->SetPrint();}

}

