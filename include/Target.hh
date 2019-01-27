#ifndef Target_H
#define Target_H 1

#include "G4Material.hh"
#include "Materials.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4ExtrudedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4AssemblyVolume.hh"
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

class Target 
{
  public:

  G4LogicalVolume *expHall_log;

  Target(Materials*);
  ~Target();
  
  G4VPhysicalVolume *Construct(G4LogicalVolume*);
  void setX(G4double);
  void setY(G4double);
  void setZ(G4double);
  void setMaterial(G4String);
  void setNStep(G4int);
  void setSourceFrame(G4String);
  void setSled(){buildSled=true;}
  void BuildSled();
  
  void Report();
  G4LogicalVolume* GetTargetLog(){return Target_log;}
  G4Box* GetTarget(){return aTarget;}
  G4VPhysicalVolume* GetTargetPlacement(){return Target_phys;}
  void setTargetReactionDepth(G4double);
  void ScaleDensity(G4double);
  //  void SetPositionZ(G4double);
  void SetPosition(G4double, G4double, G4double);
  G4double GetTargetThickness(){return Target_thickness;}
  G4ThreeVector* GetPos(){return Pos;}

private:
  // dimensions
  G4double Target_side_x;
  G4double Target_side_y;
  G4double Target_thickness;
  G4double frameInnerRadius;
  G4double frameOuterRadius;
  G4double tapeThickness;
  G4double frameSide_x;
  G4double frameSide_y;
  G4double frameThickness;

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

  Materials* materials;
  G4Material* TargetMaterial;
  G4Material* frameMaterial;
  G4Material* tapeMaterial;

  G4Material* sledMaterial;

  //default position
  G4RotationMatrix NoRot;
  G4RotationMatrix Rot;
  G4ThreeVector *Pos;

  G4Box* aTarget;
  G4Tubs* euFrame;
  G4Tubs* euTape;
  G4Tubs* csFrame;
  G4Tubs* csRing;
  G4Tubs* csTape;
  G4Box* coFrame;

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

  G4Tubs* sledFrame;
  G4Box* sledBarBoxA;
  G4Box* sledBarBoxB;
  G4Box* sledBarBoxC;
  G4AssemblyVolume* sledBar;
  G4Box* sledRunnerBox;
  G4Tubs* sledRunnerTubs;
  G4IntersectionSolid* sledRunner; 

  G4UserLimits *target_limits;

  //Number of simulation steps
  G4int NStep;

  G4String sourceFrame;
  bool buildSled;

};

#endif

