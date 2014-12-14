#ifndef Target_LH_H
#define Target_LH_H 1

#include "G4Material.hh"
#include "G4Element.hh"
#include "Materials.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4AssemblyVolume.hh"

class Target 
{
  public:

  G4LogicalVolume *expHall_log;

  Target(G4LogicalVolume*,Materials*);
  ~Target();
  
  G4VPhysicalVolume *Construct();
  void ConstructBeamTube();
  G4VPhysicalVolume *Construct200mgCell();
  G4VPhysicalVolume *Construct50mgCell();
  G4VPhysicalVolume *ConstructEmptyCell();
  G4VPhysicalVolume *ConstructNoCell();
  void ConstructWindows();
  void ConstructCryo();
  void ConstructGretaChamber();

  void setTargetCell(G4String type){ targetCellType = type; }
  G4String getTargetCell(){ return targetCellType; }
  void setBulgeThickness(G4double t){ BulgeDz = t; }
  void setWindows(){ windows = true; }
  void SetDensity(G4double d){targetDensity = d*mg/cm3;}
  void setTargetAngle(G4double a){ targetAngle = a; }
  void setX(G4double);
  void setY(G4double);
  void setZ(G4double);
  void setMaterial(G4String);
  void setNStep(G4int);
  void setSourceFrame(G4String);
  void setSled();
  void SetGreta(){ Greta = true; }
  void setCutaway(){ Cutaway = true; }
  void Report();
  G4LogicalVolume* GetTargetLog(){return Target_log;}
  G4UnionSolid* GetTarget(){return aTarget;}
  G4VPhysicalVolume* GetTargetPlacement(){return Target_phys;}
  void setTargetReactionDepth(G4double);
  void SetPositionZ(G4double);
  G4double GetTargetThickness(){return Target_thickness;}
  G4double GetDensity(){return targetDensity;}
  G4double GetTargetCellDz(){return TargetDz;}
  G4ThreeVector* GetPos(){return Pos;}

private:
  G4String targetCellType;

  G4double Target_thickness;
  G4double TargetDz;
  G4double BulgeR;
  G4double BulgeDz;

  G4bool windows;
  G4double windowThickness;

  // LH target beam tube dimensions
  G4double BeamTubeRmin;
  G4double BeamTubeRmax;
  G4double BeamTubeDz;
  G4double BeamTeeRmin;
  G4double BeamTeeRmax;
  G4double BeamTeeDz;

  //source frame dimensions
  G4double frameThickness;
  G4double frameInnerRadius;
  G4double frameOuterRadius;
  G4double tapeThickness;
  G4double tape_r;

  //target sled dimensions
  G4double sledFrameThickness;
  G4double sledFrameRmin;
  G4double sledFrameRmax;
  G4double sledRunnerThickness;
  G4double sledRunnerLength;
  G4double sledRunnerHeight;
  G4double sledRunnerRmin;
  G4double sledRunnerRmax;
  G4double sledBarThickness;
  G4double sledBarLength;
  G4double sledBarLA;
  G4double sledBarLB;
  G4double sledBarLC;
  G4double sledBarHA;
  G4double sledBarHBC;
  G4double sledBarY1;
  G4double sledBarZ1;
  G4double sledBarY2;
  G4double sledBarZ2;
  G4double sledBarY3;
  G4double sledBarZ3;
  G4double sledBarY4;
  G4double sledBarZ4;
  G4double sledBarY5;
  G4double sledBarZ5;

  // Greta Chamber dimensions
  G4double GretaChamberRmin;
  G4double GretaChamberRmax;
  G4double GretaBTrmin;
  G4double GretaBTrmax;
  G4double GretaBTDz;
  
  //elements

  //materials
  Materials* materials;
  G4Material* lH2;
  G4Material* vacuum;
  G4Material* Aluminum;
  G4Material* Copper;
  G4Material* StainlessSteel;
  G4Material* Kapton;
  G4Material* TargetMaterial;
  G4Material* frameMaterial;
  G4Material* tapeMaterial;

  G4Material* sledMaterial;

  //default position
  G4RotationMatrix NoRot;
  G4ThreeVector *Pos;
 
  //the shapes
  G4Tubs* euFrame;
  G4Tubs* euTape;
  G4Tubs* csFrame;
  G4Tubs* csRing;
  G4Tubs* csTape;
  G4Box* coFrame;

  G4Tubs* sledFrame;
  G4Box* sledBarBoxA;
  G4Box* sledBarBoxB;
  G4Box* sledBarBoxC;
  G4AssemblyVolume* sledBar;
  G4Box* sledRunnerBox;
  G4Tubs* sledRunnerTubs;
  G4IntersectionSolid* sledRunner; 

  //Assembly Volume
  G4AssemblyVolume* LHTarget;

  //Union Solid
  G4UnionSolid* aTarget;

  //logical volume
  G4LogicalVolume* Target_log;
  G4LogicalVolume* euFrame_log;
  G4LogicalVolume* euTape_log;
  G4LogicalVolume* csFrame_log;
  G4LogicalVolume* csRing_log;
  G4LogicalVolume* csTape_log;
  G4LogicalVolume* coFrame_log;

  G4LogicalVolume* sledFrame_log;
  G4LogicalVolume* sledBarBoxA_log;
  G4LogicalVolume* sledBarBoxB_log;
  G4LogicalVolume* sledBarBoxC_log;
  G4LogicalVolume* sledRunner_log;

  //physical volume
  G4VPhysicalVolume* Target_phys;
  G4VPhysicalVolume* euFrame_phys;
  G4VPhysicalVolume* euTape_phys;
  G4VPhysicalVolume* csFrame_phys;
  G4VPhysicalVolume* csRing_phys;
  G4VPhysicalVolume* csTape_phys;
  G4VPhysicalVolume* coFrame_phys;

  G4VPhysicalVolume* sledFrame_phys;
  G4VPhysicalVolume* sledRunner1_phys;
  G4VPhysicalVolume* sledRunner2_phys;

  G4UserLimits *target_limits;

  //Number of simulation steps
  G4int NStep;

  G4String sourceFrame;

  G4double targetAngle;
  G4double targetDensity;

  G4bool Greta;
  G4bool Cutaway;

};

#endif

