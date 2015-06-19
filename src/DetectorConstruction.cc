#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction()
{

  targetStatus   = false;
#ifndef SCANNING
  shellStatus    = "";
#endif
  
#ifndef LHTARGET
#ifndef SCANNING
  beamTubeStatus = false;
  gretaChamberStatus = false;
  WUChamberStatus = false;
#else
  scanningTableStatus = true;
  cloverStatus = "";
#endif
#endif

  myMessenger = new DetectorConstruction_Messenger(this);
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

  DefineMaterials();
  materials = new Materials();

  // Tracker for ions in the target

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  TrackerIon = new TrackerIonSD("IonTracker");
  TrackerIonSDMessenger = new TrackerIonSD_Messenger(TrackerIon);
  SDman->AddNewDetector( TrackerIon );

  // Tracker for gammas in GRETINA

  TrackerGamma = new TrackerGammaSD("GammaTracker");
  TrackerGammaSDMessenger = new TrackerGammaSD_Messenger(TrackerGamma);
  SDman->AddNewDetector( TrackerGamma );

  // Experimental Hall

  Experimental_Hall* ExperimentalHall = new Experimental_Hall(materials);
  ExpHall_phys = ExperimentalHall->Construct();
  ExpHall_log  = ExperimentalHall->GetLogVolume();
  ExperimentalHallMessenger = new Experimental_Hall_Messenger(ExperimentalHall);

  G4RunManager* runManager = G4RunManager::GetRunManager();
  runManager->DefineWorldVolume( ExpHall_phys );

  // Background sphere

  BackgroundSphere = new Background_Sphere(ExpHall_log,materials);
  BackgroundSphereMessenger = new Background_Sphere_Messenger(BackgroundSphere);
  BackgroundSphere->Construct();

#ifndef LHTARGET
#ifndef SCANNING
  // Beam Tube

  BeamTube = new Beam_Tube(ExpHall_log,materials);
  BeamTubeMessenger = new Beam_Tube_Messenger(BeamTube);

  // Greta Chamber

  GretaChamber = new Greta_Chamber(ExpHall_log,materials);
  GretaChamberMessenger = new Greta_Chamber_Messenger(GretaChamber);

  // WU Chamber
  WUChamber = new WU_Chamber(ExpHall_log,materials);

#else
  // Scanning Table
  scanningTable = new ScanningTable(ExpHall_log,materials);
  ScanningTableMessenger = new ScanningTable_Messenger(scanningTable);
#endif
#endif

  // Target

  aTarget = new Target(ExpHall_log,materials);
  TargetMessenger = new Target_Messenger(aTarget);

  // GRETINA

  the_Gretina_Array = new Gretina_Array();
  the_Gretina_Array_Messenger = new Gretina_Array_Messenger(the_Gretina_Array);

  return ExpHall_phys;
}

void DetectorConstruction::DefineMaterials()
{

  G4double a, z, density, temperature, pressure;
  G4double fractionmass;
  G4String name, symbol;
  G4int    nelements, natoms, nprot, nnucl;

  std::vector<G4Element*>  myElements;    // save pointers here to avoid
  std::vector<G4Material*> myMaterials;   // warnings of unused variables
  
  G4Material* Vacuum = new G4Material(name="Vacuum", z=1., a= 1.01*g/mole, density= universe_mean_density,
                                      kStateGas, temperature = 0.1*kelvin, pressure=1.0e-19*pascal);
  myMaterials.push_back(Vacuum);

  G4Element* elC  = new G4Element(name="Carbon",   symbol="C",  z=6.,  a= 12.00000*g/mole);
  myElements.push_back(elC);

  G4Element* elN  = new G4Element(name="Nitrogen", symbol="N",  z=7.,  a= 14.00674*g/mole);
  myElements.push_back(elN);

  G4Element* elO  = new G4Element(name="Oxigen",   symbol="O",  z=8.,  a= 15.9994 *g/mole);
  myElements.push_back(elO);

  G4Material* Air = new G4Material(name="Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);
  myMaterials.push_back(Air);

  G4Material* Al = new G4Material(name="Aluminium", z=13., a= 26.98154*g/mole, density= 2.70  *g/cm3);
  myMaterials.push_back(Al);

  G4Material* Si = new G4Material(name="Silicon", z=14., a= 28.0855*g/mole, density= 2.330  *g/cm3);
  myMaterials.push_back(Si);

  // Germanium isotopes
  G4Isotope* Ge70 = new G4Isotope(name="Ge70", nprot=32, nnucl=70, a=69.9242*g/mole);
  G4Isotope* Ge72 = new G4Isotope(name="Ge72", nprot=32, nnucl=72, a=71.9221*g/mole);
  G4Isotope* Ge73 = new G4Isotope(name="Ge73", nprot=32, nnucl=73, a=72.9235*g/mole);
  G4Isotope* Ge74 = new G4Isotope(name="Ge74", nprot=32, nnucl=74, a=73.9212*g/mole);
  G4Isotope* Ge76 = new G4Isotope(name="Ge76", nprot=32, nnucl=76, a=75.9214*g/mole);
  // germanium defined via its isotopes
  G4Element* elGe = new G4Element(name="Germanium",symbol="Ge", nelements=5);
  elGe->AddIsotope(Ge70, 0.2123);
  elGe->AddIsotope(Ge72, 0.2766);
  elGe->AddIsotope(Ge73, 0.0773);
  elGe->AddIsotope(Ge74, 0.3594);
  elGe->AddIsotope(Ge76, 0.0744);
  myElements.push_back(elGe);

  G4Material* Ge = new G4Material(name="Germanium", density=5.323 *g/cm3, nelements=1);
  Ge->AddElement(elGe, natoms=1);
  myMaterials.push_back(Ge);
  
  G4Element* elCu = new G4Element(name="Copper",   symbol="Cu", z=29., a=63.546    *g/mole);
  myElements.push_back(elCu);

  G4Element* elAl = new G4Element(name="Aluminum", symbol="Al", z=13., a= 26.98154*g/mole);
  myElements.push_back(elAl);

  G4Element* elMg = new G4Element(name="Magnesium",     symbol="Mg", z=12., a=24.305    *g/mole);
  myElements.push_back(elMg);

  G4Element* elMn = new G4Element(name="Manganese", symbol="Mn", z=25., a=54.938    *g/mole);
  myElements.push_back(elMn);

  G4Element* elFe = new G4Element(name="Iron",     symbol="Fe", z=26.,  a= 55.845*g/mole);
  myElements.push_back(elFe);

  G4Material* Steel = new G4Material(name="Steel", density=7.86*g/cm3, nelements=2);
  Steel->AddElement(elFe, fractionmass=98.5*perCent);
  Steel->AddElement(elC,  fractionmass= 1.5*perCent);
  myMaterials.push_back(Steel);

  // Mix of aluminum and copper behind the GRETINA crystals
  //G4Material* backWallMat = new G4Material(name="BackWallMaterial", density = 6.0*g/cm3, nelements=2);
  //  backWallMat->AddElement(elAl, 0.5);
  //  backWallMat->AddElement(elCu, 0.5);
  G4Material* backWallMat = new G4Material(name="BackWallMaterial", density = 8.96*g/cm3, nelements=1);
  backWallMat->AddElement(elCu, 1.);
  myMaterials.push_back(backWallMat);
 
  //  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

void DetectorConstruction::Placement()
{

#ifndef LHTARGET
#ifndef SCANNING
  // Beam Tube
  if( beamTubeStatus ){
    BeamTube->Construct();
  }

  // Greta Chamber
  if( gretaChamberStatus ){
    GretaChamber->Construct();
  }

  // WU Chamber
  if( WUChamberStatus ){
    WUChamber->Construct();
  }

#else
  // Scanning Table
  if( scanningTableStatus ) {
    scanningTable->Construct();
  }
  if(cloverStatus != ""){
    // Eventually fix this to place left, right, or both based on status.
    rightClover = new Clover_Detector(ExpHall_log, materials);
    rightClover->Construct();
    rightClover->MakeSensitive(TrackerGamma);
  }
#endif
#endif

  // Target
  if( targetStatus ){
    aTarget->Construct();
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
#endif
  
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
  commandName = "/ScanningTable/Construct";
  aLine = commandName.c_str();
  ScanningTableCmd = new G4UIcmdWithoutParameter(aLine, this);
  ScanningTableCmd->SetGuidance("Construct the scanning table.");
  ScanningTableCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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
#endif
  
}

DetectorConstruction_Messenger::~DetectorConstruction_Messenger()
{
  delete UpdateCmd;
  delete TargetCmd;
#ifndef SCANNING
  delete ShellCmd;
#endif
#ifndef LHTARGET
#ifndef SCANNING
  delete BeamTubeCmd;
  delete GretaChamberCmd;
  delete WUChamberCmd;
#else
  delete ScanningTableCmd;
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
  if( command == ScanningTableCmd ) {
    myTarget->SetScanningTableStatus(true);
  }
  if( command == CloverCmd ) {
    myTarget->SetCloverStatus(newValue);
  }
#endif
#endif

#ifndef SCANNING
  if( command == ShellCmd ) {
    myTarget->SetShellStatus(newValue);
  }
#endif

}
