#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "Materials.hh"
#include "Experimental_Hall.hh"
#include "Experimental_Hall_Messenger.hh"
#ifndef LHTARGET
  #include "Beam_Tube.hh"
  #include "Beam_Tube_Messenger.hh"
#else
  #include "Target_LH.hh"
#endif
#include "Target_Messenger.hh"
#include "Gretina_NSCL_Shell.hh"
#include "Gretina_Array.hh"
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

  void DefineMaterials();

  void SetTargetStatus(G4bool s){targetStatus = s;}

#ifndef LHTARGET
  void SetBeamTubeStatus(G4bool s){beamTubeStatus = s;}
#endif

  void SetShellStatus(G4String s){shellStatus = s;}

  void Placement();

private:
  DetectorConstruction_Messenger *myMessenger;

  Materials* materials;

  G4bool targetStatus;
  Target* aTarget;

#ifndef LHTARGET
  G4bool beamTubeStatus;
  Beam_Tube* BeamTube;
#endif

  G4String shellStatus;
 
  Gretina_Array*   the_Gretina_Array;

  // Logical volumes
  G4LogicalVolume* ExpHall_log;

  // Physical volumes
  G4VPhysicalVolume* ExpHall_phys;

  Experimental_Hall_Messenger* ExperimentalHallMessenger;
  Target_Messenger*    TargetMessenger;
  TrackerGammaSD* TrackerGamma;
  TrackerGammaSD_Messenger* TrackerGammaSDMessenger;
  TrackerIonSD* TrackerIon;
  TrackerIonSD_Messenger* TrackerIonSDMessenger;
  Gretina_Array_Messenger*  the_Gretina_Array_Messenger;
#ifndef LHTARGET
  Beam_Tube_Messenger* BeamTubeMessenger;
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
    G4UIcmdWithAString*          ShellCmd;
#ifndef LHTARGET
    G4UIcmdWithoutParameter*     BeamTubeCmd;
#endif
    
  public:
    void SetNewValue(G4UIcommand*, G4String);
};

#endif

