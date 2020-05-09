#ifdef SCANNING
#ifndef ScanningTable_h
#define ScanningTable_h 1

#include "G4Material.hh"
#include "Materials.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4UnitsTable.hh"

#include "CADMesh.hh"
#include "G4VisExtent.hh"

class ScanningTable 
{
  public:

  G4LogicalVolume *expHall_log;

  ScanningTable(Materials*);
  ~ScanningTable();
  
  void Construct(G4LogicalVolume*);
  void SetCADPath(G4String pth) { CADModelPath = pth; }
  void SetIncludeCloverCart() { includeCloverCart = true; }
  void SetIncludeCartFrame() { includeCartFrame = true; }
  void SetIncludeSlitMount() { includeSlitMount = true; }
  void SetIncludeSlits() { includeSlits = true; }
  void SetIncludeCollimator() { includeCollimator = true; }
  void SetIncludeCollimatorInsert() { includeCollimatorInsert = true; }
  void SetIncludeCollimatorMount() { includeCollimatorMount = true; }
  void SetIncludeShield() { includeShield = true; }
  void SetIncludeCuTarget() { includeCuTarget = true; }
  void SetControllerX(G4double value) { controllerX = value; }
  void SetControllerY(G4double value) { controllerY = value; }
  void SetControllerZ(G4double value) { controllerZ = value; }
  G4double GetControllerX() { return controllerX; }
  G4double GetControllerY() { return controllerY; }
  G4double GetControllerZ() { return controllerZ; }
  void SetCollR(G4double value)  { collimatorRadius = value; }
  void SetSlitWidth(G4double value)  { slitWidth = value; }
  void SetCloverZ(G4double value){ cloverZ = value; }
  G4double GetCloverZ() { return cloverZ + cloverOffset; }
  void Report();
  
    private:
  //flags
  G4bool includeCloverCart;
  G4bool includeCartFrame;
  G4bool includeSlits;
  G4bool includeSlitMount;
  G4bool includeCollimator;
  G4bool includeCollimatorInsert;
  G4bool includeCollimatorMount;
  G4bool includeShield;
  G4bool includeCuTarget;

  G4String CADModelPath;

  G4double controllerX;
  G4double controllerY;
  G4double controllerZ;
  G4double controllerOffsetX;
  G4double controllerOffsetY;
  G4double cloverZ;
  G4double cloverOffset;

  G4double collimatorRadius;

  G4double slitWidth;
  
  //materials
  Materials* materials;
  G4Material* material8020;
  G4Material* materialCartBase;
  G4Material* materialCartTop;
  G4Material* materialCartBack;
  G4Material* materialSlitBrackets;
  G4Material* materialSlits;
  G4Material* materialSlitAssembly;
  G4Material* materialSlitAssemblyCollimator;
  G4Material* materialSlitAssemblyPlug;
  G4Material* materialTranslation;
  G4Material* materialTranslationAssembly;
  G4Material* materialCsCollimator;
  G4Material* materialCsCollimatorMount;
  G4Material* materialClover;
  G4Material* materialCloverShield;
  G4Material* materialRollers;
  G4Material* materialCuTarget;

  //G4 Objects
  G4Box* CloverMountBase;
  G4Box* CloverMountwall;
  G4UnionSolid* CloverMount;
  G4UnionSolid* CloverMount2;
  G4Tubs* CloverMountCut;
  G4SubtractionSolid* CloverMountSub;
  G4ExtrudedSolid* CloverMountExt;
  G4ExtrudedSolid* ZSlitAssemblyTriangle;
  G4Box* CuTarget;

  //default position
  G4RotationMatrix NoRot;
  G4RotationMatrix Rot;
  G4ThreeVector *Pos0;
  G4ThreeVector *Pos1;
  G4ThreeVector *Pos2;
  G4ThreeVector *Pos3;
  G4ThreeVector *BotPos;
  G4ThreeVector *MidPos;
  G4ThreeVector *TopPos;
  //  G4ThreeVector *UsePos;
  G4ThreeVector CuTargetShift;
  G4ThreeVector CuTargetPos;
  G4ThreeVector CloverMountShift;
  G4ThreeVector CloverMount1Shift;
  G4ThreeVector ZSlitAssemblyTriangle1Shift;
  G4ThreeVector ZSlitAssemblyTriangle2Shift;
  G4RotationMatrix CloverMountRot;
  G4RotationMatrix CloverMount1Rot;
  G4RotationMatrix ZSlitAssemblyTriangleRot;
  G4ThreeVector CloverMountPos;
  G4ThreeVector CloverMount1Pos;
  G4ThreeVector ZSlitAssemblyTriangle1Pos;
  G4ThreeVector ZSlitAssemblyTriangle2Pos;

  //logical volume
  G4LogicalVolume* Cart_log;
  G4LogicalVolume* Cart8020_log;
  G4LogicalVolume* Collimator_log;
  G4LogicalVolume* CartBase_log;
  G4LogicalVolume* SlitZAssembly_log;
  G4LogicalVolume* CloverAssembly_log;
  G4LogicalVolume* CartTop_log;
  G4LogicalVolume* CuTarget_log;
  G4LogicalVolume* CartTopUBaseBottomLeft_log; 
  G4LogicalVolume* CartTopUBaseBottomMid_log;
  G4LogicalVolume* CartTopUBaseBottomRight_log;
  G4LogicalVolume* CartTopUBaseUpperLeft_log;
  G4LogicalVolume* CartTopUBaseUpperMid_log;
  G4LogicalVolume* CartTopUBaseUpperRight_log;
  G4LogicalVolume* CartTopFlange_log;
  G4LogicalVolume* CartTopPlate_log;
  G4LogicalVolume* CartTopBar_log;
  G4LogicalVolume* CartBack_log;
  G4LogicalVolume* Brake_log;
  G4LogicalVolume* ZSlit_log;
  G4LogicalVolume* ZSlitsLeftBracket_log;
  G4LogicalVolume* ZSlitsUpperLeftWall_log;
  G4LogicalVolume* ZSlitsMidLeftWall_log;
  G4LogicalVolume* ZSlitsLowerLeftWall_log;
  G4LogicalVolume* ZSlitsUpperRightWall_log;
  G4LogicalVolume* ZSlitsBottomBracket_log;
  G4LogicalVolume* ZSlitsMidRightWall_log;
  G4LogicalVolume* ZSlitsLowerRightWall_log;
  G4LogicalVolume* ZSlitsRightBracket_log;
  G4LogicalVolume* ZSlitsTopBracket_log;
  G4LogicalVolume* SlitAssembly_log;
  G4LogicalVolume* TranslateX_log;
  G4LogicalVolume* TranslateY_log;
  G4LogicalVolume* Translate_log;
  G4LogicalVolume* CollimatorBase_log;
  G4LogicalVolume* CollimatorMid1_log;
  G4LogicalVolume* CollimatorMid2_log;
  G4LogicalVolume* CollimatorMid3_log;
  G4LogicalVolume* CollimatorBody_log;
  G4LogicalVolume* CollMount_log;
  G4LogicalVolume* CloverBase_log;
  G4LogicalVolume* CloverElevator_log;
  G4LogicalVolume* CloverRight_log;
  G4LogicalVolume* CloverLeft_log;
  G4LogicalVolume* CloverAssemblyLeftBase_log;
  G4LogicalVolume* CloverAssemblyLeftLeft_log;
  G4LogicalVolume* CloverAssemblyLeftRight_log;
  G4LogicalVolume* CloverAssemblyRightBase_log;
  G4LogicalVolume* CloverAssemblyRightLeft_log;
  G4LogicalVolume* CloverAssemblyRightRight_log;
  G4LogicalVolume* CloverRightShield_log;
  G4LogicalVolume* CloverLeftShield_log;
  G4LogicalVolume* CloverMount_log;
  G4LogicalVolume* ZSlitAssemblyTriangle_log;
  
  //physical volume
  G4VPhysicalVolume* Cart_phys;
  G4VPhysicalVolume* Cart8020_phys;
  G4VPhysicalVolume* Collimator_phys;
  G4VPhysicalVolume* CartBase_phys;
  G4VPhysicalVolume* SlitZAssembly_phys;
  G4VPhysicalVolume* CloverAssembly_phys;
  G4VPhysicalVolume* CartTop_phys;
  G4VPhysicalVolume* CuTarget_phys;
  G4VPhysicalVolume* CartTopUBaseBottomLeft_phys;
  G4VPhysicalVolume* CartTopUBaseBottomMid_phys;
  G4VPhysicalVolume* CartTopUBaseBottomRight_phys;
  G4VPhysicalVolume* CartTopUBaseUpperLeft_phys;
  G4VPhysicalVolume* CartTopUBaseUpperMid_phys;
  G4VPhysicalVolume* CartTopUBaseUpperRight_phys;
  G4VPhysicalVolume* CartTopFlange_phys;
  G4VPhysicalVolume* CartTopPlate_phys;
  G4VPhysicalVolume* CartTopBar_phys;
  G4VPhysicalVolume* CartBack_phys;
  G4VPhysicalVolume* Brake_phys;
  G4VPhysicalVolume* ZSlit_phys;
  G4VPhysicalVolume* ZSlitsLeftBracket_phys;
  G4VPhysicalVolume* SlitAssembly_phys;
  G4VPhysicalVolume* ZSlitsUpperLeftWall_phys;
  G4VPhysicalVolume* ZSlitsMidLeftWall_phys;
  G4VPhysicalVolume* ZSlitsLowerLeftWall_phys;
  G4VPhysicalVolume* ZSlitsUpperRightWall_phys;
  G4VPhysicalVolume* ZSlitsBottomBracket_phys;
  G4VPhysicalVolume* ZSlitsMidRightWall_phys;
  G4VPhysicalVolume* ZSlitsLowerRightWall_phys;
  G4VPhysicalVolume* ZSlitsTopBracket_phys;
  G4VPhysicalVolume* ZSlitsRightBracket_phys;
  G4VPhysicalVolume* TranslateX_phys;
  G4VPhysicalVolume* TranslateY_phys;
  G4VPhysicalVolume* Translate_phys;
  G4VPhysicalVolume* CollimatorBase_phys;
  G4VPhysicalVolume* CollimatorMid1_phys;
  G4VPhysicalVolume* CollimatorMid2_phys;
  G4VPhysicalVolume* CollimatorMid3_phys;
  G4VPhysicalVolume* CollimatorBody_phys;
  G4VPhysicalVolume* CollMount_phys;
  G4VPhysicalVolume* CloverBase_phys;
  G4VPhysicalVolume* CloverElevator_phys;
  G4VPhysicalVolume* CloverRight_phys;
  G4VPhysicalVolume* CloverLeft_phys;
  G4VPhysicalVolume* CloverAssemblyRightBase_phys;
  G4VPhysicalVolume* CloverAssemblyRightLeft_phys;
  G4VPhysicalVolume* CloverAssemblyRightRight_phys;
  G4VPhysicalVolume* CloverAssemblyLeftBase_phys;
  G4VPhysicalVolume* CloverAssemblyLeftLeft_phys;
  G4VPhysicalVolume* CloverAssemblyLeftRight_phys;
  G4VPhysicalVolume* CloverRightShield_phys;
  G4VPhysicalVolume* CloverLeftShield_phys;
  G4VPhysicalVolume* CloverMount_phys;
  G4VPhysicalVolume* CloverMount1_phys;
  G4VPhysicalVolume* ZSlitAssemblyTriangle1_phys;
  G4VPhysicalVolume* ZSlitAssemblyTriangle2_phys;
};

#endif
#endif
