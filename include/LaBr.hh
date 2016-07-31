#ifndef LABR_H
#define LABR_H 1


#include "G4Material.hh"
#include "Materials.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

class LaBr 
{
  public:

  G4LogicalVolume *expHall_log;

  LaBr(G4LogicalVolume*,Materials*);
  ~LaBr();
  
  G4VPhysicalVolume *Construct();
  void setR(G4double);
  void setLength(G4double);
  void setMaterial(G4String);
  void Report();
  
    private:
  // dimensions
  G4double LaBr_R;
  G4double LaBr_Dz;
  G4double LaBrPMT_R;
  G4double LaBrPMT_Dz;

  //materials
  Materials* materials;
  G4Material* LaBrMaterial;
  G4Material* LaBrPMTMaterial;

  //default position
  G4RotationMatrix NoRot;
  G4ThreeVector *LaBrPos;
  G4ThreeVector *LaBrPMTPos;

  //solids
  G4Tubs* LaBr_Crys;
  G4Tubs* LaBrPMT;

  //logical volumes
  G4LogicalVolume* LaBr_log;
  G4LogicalVolume* LaBrPMT_log;
 
  //physical volumes
  G4VPhysicalVolume* LaBr_phys;
  G4VPhysicalVolume* LaBrPMT_phys;

};

#endif
