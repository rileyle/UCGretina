#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "Materials.hh"
#include "Experimental_Hall.hh"
#include "Experimental_Hall_Messenger.hh"
#include "Background_Sphere.hh"
#include "Background_Sphere_Messenger.hh"
#ifndef LHTARGET
#ifndef SCANNING
  #include "Beam_Tube.hh"
  #include "Beam_Tube_Messenger.hh"
  #include "Greta_Chamber.hh"
  #include "Greta_Chamber_Messenger.hh"
  #include "WU_Chamber.hh"
#else
  #include "ScanningTable.hh"
  #include "ScanningTable_Messenger.hh"
  #include "Clover_Detector.hh"
#endif
  #include "Target.hh"
  #include "Target_Messenger.hh"
#else
  #include "Target_LH.hh"
  #include "Target_LH_Messenger.hh"
#endif

#include "Gretina_NSCL_Shell.hh"
#include "Greta_Shell.hh"
#include "Gretina_Array.hh"
#include "S800.hh"
#include "LaBr.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "TrackerIonSD.hh"
#include "TrackerIonSD_Messenger.hh"
#include "TrackerGammaSD.hh"
#include "TrackerGammaSD_Messenger.hh"
#include "G4SDManager.hh"

class Gretina_Array;
class Gretina_Array_Messenger;
class DetectorConstruction_Messenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();
  Gretina_Array* GetGretina(){ return the_Gretina_Array;}
  TrackerGammaSD* GetGammaSD(){ return TrackerGamma;}
#ifdef LHTARGET
  G4UnionSolid* GetTarget(){return aTarget->GetTarget();}	
#else
  G4Box* GetTarget(){return aTarget->GetTarget();}
#endif
  G4VPhysicalVolume* GetTargetPlacement(){return aTarget->GetTargetPlacement();}
  G4double GetTargetThickness(){return aTarget->GetTargetThickness();}
  void setTargetReactionDepth(G4double depth){ aTarget->setTargetReactionDepth(depth);};
  G4ThreeVector* GetTargetPos(){return aTarget->GetPos();}

  G4LogicalVolume*   HallLog(){return ExpHall_log;}
  G4VPhysicalVolume* HallPhys(){return ExpHall_phys;}

  G4double GetBGSphereRmin(){return BackgroundSphere->GetRmin();}
  G4double GetBGSphereRmax(){return BackgroundSphere->GetRmax();}

  void SetTargetStatus(G4bool stat){targetStatus = stat;}
  void SetGretinaStatus(G4bool stat){gretinaStatus = stat;}

#ifndef LHTARGET
#ifndef SCANNING
  void SetBeamTubeStatus(G4bool stat){beamTubeStatus = stat;}
  void SetGretaChamberStatus(G4bool stat){gretaChamberStatus = stat;}
  void SetWUChamberStatus(G4bool stat){WUChamberStatus = stat;}
#else
  void SetCloverStatus(G4String stat){cloverStatus = stat;}
  Clover_Detector* rightClover;
  Clover_Detector* leftClover;
#endif
#endif

#ifndef SCANNING
  void SetShellStatus(G4String stat){shellStatus = stat;}
  void SetS800Status(G4bool stat){s800Status = stat;}
  void SetLaBrStatus(G4bool stat){laBrStatus = stat;}
#endif
  
  void Placement();

private:
  DetectorConstruction_Messenger *myMessenger;

  Materials* materials;

  G4bool targetStatus;
  G4bool gretinaStatus;
  Target* aTarget;

  Background_Sphere* BackgroundSphere;
  Experimental_Hall* ExperimentalHall;

#ifndef LHTARGET
#ifndef SCANNING
  G4bool beamTubeStatus;
  Beam_Tube* BeamTube;
  
  G4bool gretaChamberStatus;
  Greta_Chamber* GretaChamber;

  G4bool WUChamberStatus;
  WU_Chamber* WUChamber;
#else
  ScanningTable* scanningTable;
  G4String cloverStatus;
#endif
#endif

#ifndef SCANNING
  G4String shellStatus;
  G4bool   s800Status;
  S800*    the_S800;
  G4bool   laBrStatus;
  LaBr*    the_LaBr;
#endif
  
  Gretina_Array*   the_Gretina_Array;

  // Logical volumes
  G4LogicalVolume* ExpHall_log;

  // Physical volumes
  G4VPhysicalVolume* ExpHall_phys;

  Experimental_Hall_Messenger* ExperimentalHallMessenger;
  Target_Messenger*    TargetMessenger;
  Background_Sphere_Messenger* BackgroundSphereMessenger;
  TrackerGammaSD* TrackerGamma;
  TrackerGammaSD_Messenger* TrackerGammaSDMessenger;
  TrackerIonSD* TrackerIon;
  TrackerIonSD_Messenger* TrackerIonSDMessenger;

  Gretina_Array_Messenger*  the_Gretina_Array_Messenger;
#ifndef LHTARGET
#ifndef SCANNING
  Beam_Tube_Messenger* BeamTubeMessenger;
  Greta_Chamber_Messenger* GretaChamberMessenger;
#else
  ScanningTable_Messenger* ScanningTableMessenger;
#endif
#endif
};


#include "G4UImessenger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"

class DetectorConstruction_Messenger: public G4UImessenger
{
public:
  DetectorConstruction_Messenger(DetectorConstruction*);
  ~DetectorConstruction_Messenger();
  
private:
  DetectorConstruction*        myTarget;
  
private:
  G4UIcmdWithoutParameter*     UpdateCmd;
  G4UIcmdWithoutParameter*     TargetCmd;
  G4UIcmdWithoutParameter*     NoGretCmd;
#ifndef SCANNING
  G4UIcmdWithAString*          ShellCmd;
  G4UIcmdWithoutParameter*     S800Cmd;
  G4UIcmdWithoutParameter*     LaBrCmd;
#endif
#ifndef LHTARGET
#ifndef SCANNING
  G4UIcmdWithoutParameter*     BeamTubeCmd;
  G4UIcmdWithoutParameter*     GretaChamberCmd;
  G4UIcmdWithoutParameter*     WUChamberCmd;
#else
  G4UIcmdWithAString*          CloverCmd;
#endif
#endif
 
public:
  void SetNewValue(G4UIcommand*, G4String);
};

#endif
