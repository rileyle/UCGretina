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

  ControllerXCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetControllerX",
					    this);
  ControllerXCmd->SetGuidance("Set the X position of the scanning table controller.");
  ControllerXCmd->SetParameterName("choice", false);
  ControllerXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ControllerYCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetControllerY", this);
  ControllerYCmd->SetGuidance("Set the Y position of the scanning table controller.");
  ControllerYCmd->SetParameterName("choice", false);
  ControllerYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  ControllerZCmd = new G4UIcmdWithADoubleAndUnit("/ScanningTable/SetControllerZ", this);
  ControllerZCmd->SetGuidance("Set the vertical position of the scanning table slits.");
  ControllerZCmd->SetParameterName("choice", false);
  ControllerZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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
  delete ControllerXCmd;
  delete ControllerYCmd;
  delete ControllerZCmd;
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
  if ( command == ControllerXCmd )
    { scanningTable->SetControllerX(ControllerXCmd->GetNewDoubleValue(newValue)); }
  if ( command == ControllerYCmd )
    { scanningTable->SetControllerY(ControllerYCmd->GetNewDoubleValue(newValue)); }
  if ( command == ControllerZCmd )
    { scanningTable->SetControllerZ(ControllerZCmd->GetNewDoubleValue(newValue)); }
  if ( command == CloverZCmd )
    { scanningTable->SetCloverZ(CloverZCmd->GetNewDoubleValue(newValue)); }
  if ( command == CollRCmd )
    { scanningTable->SetCollR(CollRCmd->GetNewDoubleValue(newValue)); }
  if ( command == SlitWCmd )
    { scanningTable->SetSlitWidth(SlitWCmd->GetNewDoubleValue(newValue)); }
}
#endif
