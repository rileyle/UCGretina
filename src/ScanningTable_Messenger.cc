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
  CloverCartCmd->SetGuidance("Include clover cart.");

  CartFrameCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeCartFrame", this);
  CartFrameCmd->SetGuidance("Include the cart frame and GRETINA module mount.");

  SlitMountCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeSlitMount", this);
  SlitMountCmd->SetGuidance("Include the slit assembly mount.");

  CollimatorCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeCollimator", this);
  CollimatorCmd->SetGuidance("Include the collimator.");

  CollimatorInsertCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeCollimatorInsert", this);
  CollimatorInsertCmd->SetGuidance("Include the collimator insert.");

  CollimatorMountCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeCollimatorMount", this);
  CollimatorMountCmd->SetGuidance("Include the collimator mount.");

  ShieldCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeShields", this);
  ShieldCmd->SetGuidance("Include the BGO anti-Compton shields.");

  CuTargetCmd = new G4UIcmdWithoutParameter("/ScanningTable/IncludeCuTarget", this);
  CuTargetCmd->SetGuidance("Include the Cu Target.");

  XShiftCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetXShift", this);
  XShiftCmd->SetGuidance("Set the horizontal (X) shift of the source from nominal.");
  XShiftCmd->SetParameterName("choice", false);
  XShiftCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  YShiftCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetYShift", this);
  YShiftCmd->SetGuidance("Set the horizontal (Y) shift of the source from nominal.");
  YShiftCmd->SetParameterName("choice", false);
  YShiftCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ZShiftCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetZShift", this);
  ZShiftCmd->SetGuidance("Set the vertical shift of the slits from nominal.");
  ZShiftCmd->SetParameterName("choice", false);
  ZShiftCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CloverZCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetCloverZ", this);
  CloverZCmd->SetGuidance("Set the vertical position of the Clover(s).");
  CloverZCmd->SetParameterName("choice", false);
  CloverZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CollRCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetCollimatorRadius", this);
  CollRCmd->SetGuidance("Set the inner radius of the collimator insert.");
  CollRCmd->SetParameterName("choice", false);
  CollRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SlitWCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetSlitWidth", this);
  SlitWCmd->SetGuidance("Set the slit width.");
  SlitWCmd->SetParameterName("choice", false);
  SlitWCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

ScanningTable_Messenger::~ScanningTable_Messenger()
{
  delete ScanningTableDir;
  delete CloverCartCmd;
  delete CADPathCmd;
  delete CartFrameCmd;
  delete SlitMountCmd;
  delete CollimatorCmd;
  delete CollimatorInsertCmd;
  delete CollimatorMountCmd;
  delete ShieldCmd;
  delete CuTargetCmd;
  delete XShiftCmd;
  delete YShiftCmd;
  delete ZShiftCmd;
  delete CloverZCmd;
  delete CollRCmd;
  delete SlitWCmd;
}

void ScanningTable_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if ( command == CloverCartCmd )
    { scanningTable->SetIncludeCloverCart(); }
  if ( command == CADPathCmd )
    { scanningTable->SetCADPath(newValue); }
  if ( command == CartFrameCmd )
    { scanningTable->SetIncludeCartFrame(); }
  if ( command == SlitMountCmd )
    { scanningTable->SetIncludeSlitMount(); }
  if ( command == CollimatorCmd )
    { scanningTable->SetIncludeCollimator(); }
  if ( command == CollimatorInsertCmd )
    { scanningTable->SetIncludeCollimatorInsert(); }
  if ( command == CollimatorMountCmd )
    { scanningTable->SetIncludeCollimatorMount(); }
  if ( command == ShieldCmd )
    { scanningTable->SetIncludeShield(); }
  if ( command == CuTargetCmd )
    { scanningTable->SetIncludeCuTarget(); }
  if ( command == XShiftCmd )
    { scanningTable->SetXShift(XShiftCmd->GetNewDoubleValue(newValue)); }
  if ( command == YShiftCmd )
    { scanningTable->SetYShift(YShiftCmd->GetNewDoubleValue(newValue)); }
  if ( command == ZShiftCmd )
    { scanningTable->SetZShift(ZShiftCmd->GetNewDoubleValue(newValue)); }
  if ( command == CloverZCmd )
    { scanningTable->SetCloverZ(CloverZCmd->GetNewDoubleValue(newValue)); }
  if ( command == CollRCmd )
    { scanningTable->SetCollR(CollRCmd->GetNewDoubleValue(newValue)); }
  if ( command == SlitWCmd )
    { scanningTable->SetSlitWidth(SlitWCmd->GetNewDoubleValue(newValue)); }
}
#endif
