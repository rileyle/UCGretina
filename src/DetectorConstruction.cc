#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction()
{

  targetStatus   = false;
#ifndef SCANNING
  shellStatus    = "";
  s800Status     = false;
  laBrStatus     = false;
#endif
  
#ifndef LHTARGET
#ifndef SCANNING
  beamTubeStatus = false;
  gretaChamberStatus = false;
  WUChamberStatus = false;
#else
  cloverStatus = "";
#endif
#endif

  gretinaStatus = true;

  myMessenger = new DetectorConstruction_Messenger(this);

  //  DefineMaterials();
  materials = new Materials();

  TrackerIon = new TrackerIonSD("IonTracker");
  TrackerIonSDMessenger = new TrackerIonSD_Messenger(TrackerIon);

  TrackerGamma = new TrackerGammaSD("GammaTracker");
  TrackerGammaSDMessenger = new TrackerGammaSD_Messenger(TrackerGamma);

  ExperimentalHall = new Experimental_Hall(materials);
  ExperimentalHallMessenger = new Experimental_Hall_Messenger(ExperimentalHall);
  
  // Background sphere
  BackgroundSphere = new Background_Sphere(materials);
  BackgroundSphereMessenger = new Background_Sphere_Messenger(BackgroundSphere);

#ifndef LHTARGET
#ifndef SCANNING
  // Beam Tube

  BeamTube = new Beam_Tube(materials);
  BeamTubeMessenger = new Beam_Tube_Messenger(BeamTube);

  // Greta Chamber

  GretaChamber = new Greta_Chamber(materials);
  GretaChamberMessenger = new Greta_Chamber_Messenger(GretaChamber);

  // WU Chamber
  WUChamber = new WU_Chamber(materials);

#else
  // Scanning Table
  scanningTable = new ScanningTable(materials);
  ScanningTableMessenger = new ScanningTable_Messenger(scanningTable);
#endif
#endif

#ifndef SCANNING  
  // S800 Quadrupole
  the_S800 = new S800(materials);
  the_LaBr = new LaBr(materials);
#endif
  
  // Target

  aTarget = new Target(materials);
  TargetMessenger = new Target_Messenger(aTarget);

  // GRETINA

  the_Gretina_Array = new Gretina_Array();
  the_Gretina_Array_Messenger = new Gretina_Array_Messenger(the_Gretina_Array);

}

DetectorConstruction::~DetectorConstruction()
{
  delete ExperimentalHallMessenger;
  delete TargetMessenger;
  delete TrackerIonSDMessenger;
  delete TrackerGammaSDMessenger;
#ifndef LHTARGET
#ifndef SCANNING
  delete BeamTubeMessenger;
  delete GretaChamberMessenger;
#else
  delete ScanningTableMessenger;
#endif
#endif
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  // Tracker for ions in the target

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SDman->AddNewDetector( TrackerIon );

  // Tracker for gammas in GRETINA

  SDman->AddNewDetector( TrackerGamma );

  // Experimental Hall

  ExpHall_phys = ExperimentalHall->Construct();
  ExpHall_log  = ExperimentalHall->GetLogVolume();

  G4RunManager* runManager = G4RunManager::GetRunManager();
  runManager->DefineWorldVolume( ExpHall_phys );

  BackgroundSphere->Construct(ExpHall_log);

  Placement();

  return ExpHall_phys;
}

void DetectorConstruction::DefineMaterials()
{
  // Remove this method. All materials are now defined by the Materials class.
}

void DetectorConstruction::Placement()
{

#ifndef LHTARGET
#ifndef SCANNING
  // Beam Tube
  if( beamTubeStatus ){
    BeamTube->Construct(ExpHall_log);
  }

  // Greta Chamber
  if( gretaChamberStatus ){
    GretaChamber->Construct(ExpHall_log);
  }

  // WU Chamber
  if( WUChamberStatus ){
    WUChamber->Construct(ExpHall_log);
  }

#else
  // Scanning Table
  scanningTable->Construct(ExpHall_log);
  if( cloverStatus == "left"  || cloverStatus == "both" ){
    G4cout << "Constructing left clover detector." << G4endl;
    leftClover = new Clover_Detector(ExpHall_log, materials, "left");
    leftClover->setY(scanningTable->GetCloverZ());
    leftClover->Construct();
    leftClover->MakeSensitive(TrackerGamma);
  }
  if( cloverStatus == "right" || cloverStatus == "both" ){
    G4cout << "Constructing right clover detector." << G4endl;
    rightClover = new Clover_Detector(ExpHall_log, materials, "right");
    rightClover->setY(scanningTable->GetCloverZ());
    rightClover->Construct();
    rightClover->MakeSensitive(TrackerGamma);
  }
#endif
#endif

  // Target
  if( targetStatus ){
    aTarget->Construct(ExpHall_log);
    aTarget->GetTargetLog()->SetSensitiveDetector(TrackerIon);
  }

#ifndef SCANNING
  if( shellStatus == "full"  ||
      shellStatus == "north" ||
      shellStatus == "south"){
    Gretina_NSCL_Shell* Shell = new Gretina_NSCL_Shell();
    Shell->Placement(shellStatus);
  } else if ( shellStatus == "Greta" ||
	      shellStatus == "GretaLH" || 
	      shellStatus == "Greta_North" ||
	      shellStatus == "Greta_South" ||
	      shellStatus == "GretaLH_North" ||
	      shellStatus == "GretaLH_South" ){
    Greta_Shell* Shell = new Greta_Shell();
    Shell->Placement(shellStatus);
  }
  if( s800Status ){
    the_S800->Construct(ExpHall_log);
  }
  if( laBrStatus ){
    the_LaBr->Construct(ExpHall_log);
  }
#endif

  if(gretinaStatus)
    the_Gretina_Array->Placement();

}

///////////////////
// The Messenger
///////////////////

DetectorConstruction_Messenger::DetectorConstruction_Messenger(DetectorConstruction* pTarget):myTarget(pTarget)
{

  const char *aLine;
  G4String commandName;

  commandName = "/Gretina/update";
  aLine = commandName.c_str();
  UpdateCmd = new G4UIcmdWithoutParameter(aLine, this);
  UpdateCmd->SetGuidance("Update geometry.");
  UpdateCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/Target/Construct";
  aLine = commandName.c_str();
  TargetCmd = new G4UIcmdWithoutParameter(aLine, this);
  TargetCmd->SetGuidance("Construct the target.");
  TargetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/Gretina/NoDetectors";
  aLine = commandName.c_str();
  NoGretCmd = new G4UIcmdWithoutParameter(aLine, this);
  NoGretCmd->SetGuidance("Omit the GRETINA detectors.");
  NoGretCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

#ifndef LHTARGET
#ifndef SCANNING
  commandName = "/BeamTube/Construct";
  aLine = commandName.c_str();
  BeamTubeCmd = new G4UIcmdWithoutParameter(aLine, this);
  BeamTubeCmd->SetGuidance("Construct the beam tube.");
  BeamTubeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/GretaChamber/Construct";
  aLine = commandName.c_str();
  GretaChamberCmd = new G4UIcmdWithoutParameter(aLine, this);
  GretaChamberCmd->SetGuidance("Construct the Greta chamber.");
  GretaChamberCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/WUChamber/Construct";
  aLine = commandName.c_str();
  WUChamberCmd = new G4UIcmdWithoutParameter(aLine, this);
  WUChamberCmd->SetGuidance("Construct the WU chamber.");
  WUChamberCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

#else
  commandName = "/ScanningTable/Clover";
  aLine = commandName.c_str();
  CloverCmd = new G4UIcmdWithAString(aLine, this);
  CloverCmd->SetGuidance("Construct the clover detector(s) (left/right/both).");
  CloverCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
#endif
#endif

#ifndef SCANNING
  commandName = "/Gretina/Shell";
  aLine = commandName.c_str();
  ShellCmd = new G4UIcmdWithAString(aLine, this);
  ShellCmd->SetGuidance("Construct the GRETINA mounting shell (full/north/south).");
  ShellCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/Gretina/S800";
  aLine = commandName.c_str();
  S800Cmd = new G4UIcmdWithoutParameter(aLine, this);
  S800Cmd->SetGuidance("Construct the S800 Quadrupole.");
  S800Cmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  commandName = "/Gretina/LaBr";
  aLine = commandName.c_str();
  LaBrCmd = new G4UIcmdWithoutParameter(aLine, this);
  LaBrCmd->SetGuidance("Construct the LaBr.");
  LaBrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
#endif
  
}

DetectorConstruction_Messenger::~DetectorConstruction_Messenger()
{
  delete UpdateCmd;
  delete TargetCmd;
  delete NoGretCmd;
#ifndef SCANNING
  delete ShellCmd;
  delete S800Cmd;
  delete LaBrCmd;
#endif
#ifndef LHTARGET
#ifndef SCANNING
  delete BeamTubeCmd;
  delete GretaChamberCmd;
  delete WUChamberCmd;
#else
  delete CloverCmd;
#endif
#endif
}

void DetectorConstruction_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == UpdateCmd ) {
    myTarget->Placement();
  } 
  if( command == TargetCmd ) {
    myTarget->SetTargetStatus(true);
  } 
  if( command == NoGretCmd ) {
    myTarget->SetGretinaStatus(false);
  }
#ifndef LHTARGET
#ifndef SCANNING
  if( command == BeamTubeCmd ) {
    myTarget->SetBeamTubeStatus(true);
  } 
  if( command == GretaChamberCmd ) {
    myTarget->SetGretaChamberStatus(true);
  } 
  if( command == WUChamberCmd ) {
    myTarget->SetWUChamberStatus(true);
  }
#else
  if( command == CloverCmd ) {
    myTarget->SetCloverStatus(newValue);
  }
#endif
#endif

#ifndef SCANNING
  if( command == ShellCmd ) {
    myTarget->SetShellStatus(newValue);
  }
  if( command == S800Cmd ) {
    myTarget->SetS800Status(true);
  } 
  if( command == LaBrCmd ) {
    myTarget->SetLaBrStatus(true);
  } 

#endif

}
