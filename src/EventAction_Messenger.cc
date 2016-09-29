#include "EventAction_Messenger.hh"


EventAction_Messenger::EventAction_Messenger(EventAction* EA):theEventAction(EA) 
{ 

  OutputDir = new G4UIdirectory("/Output/");
  OutputDir->SetGuidance("Event by event output control.");

  OutFileCmd = new G4UIcmdWithAString("/Output/Filename",this);
  OutFileCmd->SetGuidance("Output file name.");

  OutDetCmd = new G4UIcmdWithoutParameter("/Output/DetectorsOnly",this);
  OutDetCmd->SetGuidance("Only print detected gamma-ray information.");

  PosResCmd = new G4UIcmdWithADoubleAndUnit("/Gretina/detector/PositionResolution", this);
  PosResCmd->SetGuidance("Set position resolution (default = 0).");

  ThreshECmd = new G4UIcmdWithADoubleAndUnit("/Gretina/detector/ThresholdEnergy", this);
  ThreshECmd->SetGuidance("Set crystal low-energy threshold (default = 0).");

  ThreshDECmd = new G4UIcmdWithADoubleAndUnit("/Gretina/detector/ThresholdDE", this);
  ThreshDECmd->SetGuidance("Set crystal low energy threshold width paramter (default = 0.001 keV).");
  
  Mode2Dir = new G4UIdirectory("/Mode2/");
  Mode2Dir->SetGuidance("Parameters for simulating Mode2 data.");

  Mode2FileCmd = new G4UIcmdWithAString("/Mode2/Filename",this);
  Mode2FileCmd->SetGuidance("Mode 2 output file name");

  crmatCmd = new G4UIcmdWithAString("/Mode2/crmatFile",this);
  crmatCmd->SetGuidance("Use the crystal-frame to world-frame transformations in the specified file for Mode2 data (expected in crystal coordinates).");

  crysCmd = new G4UIcmdWithoutParameter("/Mode2/crystalXforms",this);
  crysCmd->SetGuidance("Use internal transformations from world coordinates to crystal frames for Mode2 data (expected in crystal coordinates).");

  coordsCmd = new G4UIcmdWithoutParameter("/Mode2/GretinaCoords",this);
  coordsCmd->SetGuidance("Write interaction points in GRETINA coordinate system --- x down, z beam (standard for Mode2 data).");

  PackResCmd = new G4UIcmdWithADoubleAndUnit("/Mode2/PackingRes",this);
  PackResCmd->SetGuidance("Packing resolution for simulating Mode2 data.");

  S800KECmd = new G4UIcmdWithADoubleAndUnit("/Mode2/S800KE",this);
  S800KECmd->SetGuidance("Average KE for computing dta for writing S800 events.");

  printCmd = new G4UIcmdWithoutParameter("/Mode2/Print",this);
  printCmd->SetGuidance("Write mode 2 event information to stdout.");

}



EventAction_Messenger::~EventAction_Messenger()

{  

  delete OutputDir;
  delete OutFileCmd;
  delete OutDetCmd;
  delete PosResCmd;
  delete ThreshECmd;
  delete ThreshDECmd;
  delete Mode2Dir;
  delete Mode2FileCmd;
  delete crmatCmd;
  delete crysCmd;
  delete coordsCmd;
  delete PackResCmd;
  delete S800KECmd;
  delete printCmd;

}



void EventAction_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == OutFileCmd )
    {theEventAction->SetOutFile(newValue);}
  if( command == OutDetCmd )
    {theEventAction->SetOutDetsOnly();}
  if( command == PosResCmd )
    {theEventAction->SetPosRes(PosResCmd->GetNewDoubleValue(newValue));}
  if( command == ThreshECmd )
    {theEventAction->SetThreshE(ThreshECmd->GetNewDoubleValue(newValue));}
  if( command == ThreshDECmd )
    {theEventAction->SetThreshDE(ThreshDECmd->GetNewDoubleValue(newValue));}
  if( command == Mode2FileCmd )
    {theEventAction->SetMode2File(newValue);}
  if( command == crmatCmd )
    {theEventAction->SetCrmatFile(newValue);}
  if( command == crysCmd )
    {theEventAction->SetCrystalXforms();}
  if( command == coordsCmd )
    {theEventAction->SetGretinaCoords();}
  if( command == PackResCmd )
    {theEventAction->SetPackRes(PackResCmd->GetNewDoubleValue(newValue));}
  if( command == S800KECmd )
    {theEventAction->SetS800KE(S800KECmd->GetNewDoubleValue(newValue));}
  if( command == printCmd )
    {theEventAction->SetPrint();}

}

