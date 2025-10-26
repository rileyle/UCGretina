#include "EventAction_Messenger.hh"


EventAction_Messenger::EventAction_Messenger(EventAction* EA):theEventAction(EA) 
{ 

  OutputDir = new G4UIdirectory("/Output/");
  OutputDir->SetGuidance("Event by event output control.");

  OutFileCmd = new G4UIcmdWithAString("/Output/Filename",this);
  OutFileCmd->SetGuidance("Output file name.");

  OutDetCmd = new G4UIcmdWithoutParameter("/Output/DetectorsOnly",this);
  OutDetCmd->SetGuidance("Only print detected gamma-ray information.");

  ThreshECmd = new G4UIcmdWithADoubleAndUnit("/Gretina/detector/ThresholdEnergy", this);
  ThreshECmd->SetGuidance("Set crystal low-energy threshold (default = 0).");

  ThreshDECmd = new G4UIcmdWithADoubleAndUnit("/Gretina/detector/ThresholdDE", this);
  ThreshDECmd->SetGuidance("Set crystal low energy threshold width paramter (default = 0.001 keV).");
  
  Mode2Dir = new G4UIdirectory("/Mode2/");
  Mode2Dir->SetGuidance("Parameters for simulating Mode2 data.");

  Mode2FileCmd = new G4UIcmdWithAString("/Mode2/Filename",this);
  Mode2FileCmd->SetGuidance("Mode 2 output file name");

  cacheDir = new G4UIdirectory("/Cache/");
  cacheOutputFileCmd = new G4UIcmdWithAString("/Cache/Output",this);
  cacheOutputFileCmd -> SetGuidance("sets the EventAction Cache Ouput file");
  
  cacheInputFileCmd = new G4UIcmdWithAString("/Cache/Input",this);
  cacheInputFileCmd -> SetGuidance("sets the EventAction Cache Input file");
  
  cacheHalfLifeCmd = new G4UIcmdWithADoubleAndUnit("/Cache/HalfLife",this);
  cacheHalfLifeCmd -> SetGuidance("sets the Half-Life for the particles in the cache simulation");
  cacheGammaEnergyCmd = new G4UIcmdWithADoubleAndUnit("/Cache/GammaEnergy",this);
  cacheGammaEnergyCmd -> SetGuidance("sets the Gamma Energy for the particles in the cache simulation");
  cacheZOffsetCmd = new G4UIcmdWithADoubleAndUnit("/Cache/ZOffset",this);
  cacheZOffsetCmd -> SetGuidance("sets the Z-Offset for the target in the cache simulation");

  CacheAngDistDir = new G4UIdirectory("/Cache/AngularDistribution/");
  CacheAngDistA0Cmd = new G4UIcmdWithADouble("/Cache/AngularDistribution/a0",this);
  CacheAngDistA0Cmd->SetGuidance("a0 parameter of the gamma-ray angular distribution for cached events.");

  CacheAngDistA2Cmd = new G4UIcmdWithADouble("/Cache/AngularDistribution/a2",this);
  CacheAngDistA2Cmd->SetGuidance("a2 parameter of the gamma-ray angular distribution for cached events.");

  CacheAngDistA4Cmd = new G4UIcmdWithADouble("/Cache/AngularDistribution/a4",this);
  CacheAngDistA4Cmd->SetGuidance("a4 parameter of the gamma-ray angular distribution for cached events.");

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

  allS800Cmd = new G4UIcmdWithoutParameter("/Mode2/AllS800",this);
  allS800Cmd->SetGuidance("Write S800 data for all events.");

  TimeSortCmd = new G4UIcmdWithoutParameter("/Mode2/TimeSort",this);
  TimeSortCmd->SetGuidance("Sort interaction points in each crystal by time.");

  printCmd = new G4UIcmdWithoutParameter("/Mode2/Print",this);
  printCmd->SetGuidance("Write mode 2 event information to stdout.");

}



EventAction_Messenger::~EventAction_Messenger()

{  

  delete OutputDir;
  delete OutFileCmd;
  delete OutDetCmd;
  delete ThreshECmd;
  delete ThreshDECmd;
  delete Mode2Dir;
  delete Mode2FileCmd;
  delete cacheOutputFileCmd;
  delete cacheInputFileCmd;
  delete cacheHalfLifeCmd;
  delete cacheGammaEnergyCmd;
  delete cacheZOffsetCmd;
  delete CacheAngDistDir;
  delete CacheAngDistA0Cmd;
  delete CacheAngDistA2Cmd;
  delete CacheAngDistA4Cmd;
  delete crmatCmd;
  delete crysCmd;
  delete coordsCmd;
  delete PackResCmd;
  delete S800KECmd;
  delete allS800Cmd;
  delete TimeSortCmd;
  delete printCmd;

}



void EventAction_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == OutFileCmd )
    {theEventAction->SetOutFile(newValue);}
  if( command == OutDetCmd )
    {theEventAction->SetOutDetsOnly();}
  if( command == ThreshECmd )
    {theEventAction->SetThreshE(ThreshECmd->GetNewDoubleValue(newValue));}
  if( command == ThreshDECmd )
    {theEventAction->SetThreshDE(ThreshDECmd->GetNewDoubleValue(newValue));}
  if( command == Mode2FileCmd )
    {theEventAction->SetMode2File(newValue);}
  if( command ==  cacheOutputFileCmd )
    {theEventAction->openCacheOutputFile(newValue);}
  if( command == cacheInputFileCmd )
    {theEventAction->openCacheInputFile(newValue);}
  if( command == cacheHalfLifeCmd )
    {theEventAction->SetCacheHalfLife(cacheHalfLifeCmd->GetNewDoubleValue(newValue));}
  if( command == cacheGammaEnergyCmd )
    {theEventAction->SetCacheGammaEnergy(cacheGammaEnergyCmd->GetNewDoubleValue(newValue));}
  if( command == cacheZOffsetCmd )
    {theEventAction->SetCacheZOffset( cacheZOffsetCmd->GetNewDoubleValue(newValue));}
  if( command == CacheAngDistA0Cmd )
    {theEventAction->SetCacheAngDistA0(CacheAngDistA0Cmd->GetNewDoubleValue(newValue));}
  if( command == CacheAngDistA2Cmd )
    {theEventAction->SetCacheAngDistA2(CacheAngDistA2Cmd->GetNewDoubleValue(newValue));}
  if( command == CacheAngDistA4Cmd )
    {theEventAction->SetCacheAngDistA4(CacheAngDistA4Cmd->GetNewDoubleValue(newValue));}
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
  if( command == allS800Cmd )
    {theEventAction->SetAllS800();}
  if( command == TimeSortCmd )
    {theEventAction->SetTimeSort();}
  if( command == printCmd )
    {theEventAction->SetPrint();}

}

