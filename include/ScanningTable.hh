#ifdef SCANNING
#ifndef ScanningTable_H
#define ScanningTable_H 1

#include "G4Material.hh"
#include "Materials.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
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

  ScanningTable(G4LogicalVolume*,Materials*);
  ~ScanningTable();
  
  G4VPhysicalVolume *Construct();
  void SetCADPath(G4String pth) { CADModelPath = pth; }
  void SetIncludeCloverCart() { includeCloverCart = true; }
  void SetIncludeCartFrame() { includeCartFrame = true; }
  void SetIncludeShield() { includeShield = true; }
  void SetXShift(G4double value) { xShift = value; }
  void SetYShift(G4double value) { yShift = value; }
  void SetZShift(G4double value) { zShift = value; }
  void Report();
  
    private:
  //flags
  G4bool includeCloverCart;
  G4bool includeCartFrame;
  G4bool includeShield;

  G4String CADModelPath;

  G4double xShift;
  G4double yShift;
  G4double zShift;
  
  //materials
  Materials* materials;
  G4Material* material8020;
  G4Material* materialCartBase;
  G4Material* materialCartTop;
  G4Material* materialCartBack;
  G4Material* materialSlits;
  G4Material* materialSlitAssembly;
  G4Material* materialTranslation;
  G4Material* materialTranslationAssembly;
  G4Material* materialCsCollimator;
  G4Material* materialClover;
  G4Material* materialCloverShield;

  //default position
  G4RotationMatrix NoRot;
  G4RotationMatrix Rot;
  G4ThreeVector *Pos0;

  //logical volume
  G4LogicalVolume* Cart8020_log;
  G4LogicalVolume* CartBase_log;
  G4LogicalVolume* CartTop_log;
  G4LogicalVolume* CartBack_log;
  G4LogicalVolume* Brake_log;
  G4LogicalVolume* Slits_log;
  G4LogicalVolume* SlitAssembly_log;
  G4LogicalVolume* TranslateX_log;
  G4LogicalVolume* TranslateY_log;
  G4LogicalVolume* Translate_log;
  G4LogicalVolume* CsCollimatorBase_log;
  G4LogicalVolume* CsCollimatorBody_log;
  G4LogicalVolume* CsCollimatorPlug_log;
  G4LogicalVolume* CloverBase_log;
  G4LogicalVolume* CloverElevator_log;
  G4LogicalVolume* CloverRight_log;
  G4LogicalVolume* CloverLeft_log;
  G4LogicalVolume* CloverAssemblyLeft_log;
  G4LogicalVolume* CloverAssemblyRight_log;
  G4LogicalVolume* CloverRightShield_log;
  G4LogicalVolume* CloverLeftShield_log;
  
  //physical volume
  G4VPhysicalVolume* Cart8020_phys;
  G4VPhysicalVolume* CartBase_phys;
  G4VPhysicalVolume* CartTop_phys;
  G4VPhysicalVolume* CartBack_phys;
  G4VPhysicalVolume* Brake_phys;
  G4VPhysicalVolume* Slits_phys;
  G4VPhysicalVolume* SlitAssembly_phys;
  G4VPhysicalVolume* TranslateX_phys;
  G4VPhysicalVolume* TranslateY_phys;
  G4VPhysicalVolume* Translate_phys;
  G4VPhysicalVolume* CsCollimatorBase_phys;
  G4VPhysicalVolume* CsCollimatorBody_phys;
  G4VPhysicalVolume* CsCollimatorPlug_phys;
  G4VPhysicalVolume* CloverBase_phys;
  G4VPhysicalVolume* CloverElevator_phys;
  G4VPhysicalVolume* CloverRight_phys;
  G4VPhysicalVolume* CloverLeft_phys;
  G4VPhysicalVolume* CloverAssemblyRight_phys;
  G4VPhysicalVolume* CloverAssemblyLeft_phys;
  G4VPhysicalVolume* CloverRightShield_phys;
  G4VPhysicalVolume* CloverLeftShield_phys;

};

#endif
#endif
