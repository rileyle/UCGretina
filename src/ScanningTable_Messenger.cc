#ifdef SCANNING
#include "ScanningTable_Messenger.hh"

ScanningTable_Messenger::ScanningTable_Messenger(ScanningTable* ST)
:scanningTable(ST)
{ 
  ScanningTableDir = new G4UIdirectory("/ScanningTable/");
  ScanningTableDir->SetGuidance("Scanning table control.");

  CADPathCmd = new G4UIcmdWithAString("/ScanningTable/CADModelPath", this);
  CADPathCmd->SetGuidance("Path to STL files.");

  CloverCartCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeCloverCart", this);
  CloverCartCmd->SetGuidance("Include clover secondary cart construction.");

  CartFrameCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeCartFrame", this);
  CartFrameCmd->SetGuidance("Include cart frame construction.");

  ShieldCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeShields", this);
  ShieldCmd->SetGuidance("Include anti-Compton shields.");

  XShiftCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetXShift", this);
  XShiftCmd->SetGuidance("Set the horizontal (X) shift from nominal.");
  XShiftCmd->SetParameterName("choice", false);
  XShiftCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  YShiftCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetYShift", this);
  YShiftCmd->SetGuidance("Set the horizontal (Y) shift from nominal.");
  YShiftCmd->SetParameterName("choice", false);
  YShiftCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ZShiftCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetZShift", this);
  ZShiftCmd->SetGuidance("Set the vertical shift from nominal.");
  ZShiftCmd->SetParameterName("choice", false);
  ZShiftCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

ScanningTable_Messenger::~ScanningTable_Messenger()
{
  delete ScanningTableDir;
  delete CloverCartCmd;
  delete CADPathCmd;
  delete CartFrameCmd;
  delete ShieldCmd;
  delete XShiftCmd;
  delete YShiftCmd;
  delete ZShiftCmd;
}

void ScanningTable_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == CloverCartCmd )
    { scanningTable->SetIncludeCloverCart(); }
  if ( command == CADPathCmd )
    { scanningTable->SetCADPath(newValue); }
  if ( command == CartFrameCmd )
    { scanningTable->SetIncludeCartFrame(); }
  if ( command == ShieldCmd )
    { scanningTable->SetIncludeShield(); }
  if ( command == XShiftCmd )
    { scanningTable->SetXShift(XShiftCmd->GetNewDoubleValue(newValue)); }
  if ( command == YShiftCmd )
    { scanningTable->SetYShift(YShiftCmd->GetNewDoubleValue(newValue)); }
  if ( command == ZShiftCmd )
    { scanningTable->SetZShift(ZShiftCmd->GetNewDoubleValue(newValue)); }
  
}
#endif
