#ifndef S800_H
#define S800_H 1


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

class S800 
{
  public:

  G4LogicalVolume *expHall_log;

  S800(G4LogicalVolume*,Materials*);
  ~S800();
  
  G4VPhysicalVolume *Construct();
  void setR(G4double);
  void setLength(G4double);
  void setMaterial(G4String);
  void Report();
  
    private:
  // dimensions
  G4double S800_R;
  G4double S800_Dz;

  //materials
  Materials* materials;
  G4Material* S800Material;

  //default position
  G4RotationMatrix NoRot;
  G4ThreeVector *Pos0;

  //the tube
  G4Tubs* S800_Quad;

  //logical volume
  G4LogicalVolume* S800_log;
 
  //physical volume
  G4VPhysicalVolume* S800_phys;

};

#endif
