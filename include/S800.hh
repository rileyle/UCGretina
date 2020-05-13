#ifndef S800_H
#define S800_H 1


#include "G4Material.hh"
#include "Materials.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
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

  S800();
  ~S800();
  
  G4VPhysicalVolume *Construct(G4LogicalVolume*);
  void setR(G4double);
  void setLength(G4double);
  void setMaterial(G4String);
  void Report();
  
    private:
  // dimensions
  G4double S800_R;
  G4double S800_Dz;
  G4double S800GateValve_Dx;
  G4double S800GateValve_Dy;
  G4double S800GateValve_Dz;
  G4double S800GateValve_yOffset;
  G4double S800GateValve_zOffset;
  G4double S800_Offset;

  //materials
  G4Material* S800Material;

  //default position
  G4RotationMatrix NoRot;
  G4ThreeVector *S800_Pos;
  G4ThreeVector *S800GateValve_Pos;

  // solids
  G4Tubs*             S800_Quad;
  G4Box*              S800GateValve_Box;
  G4Tubs*             S800GateValve_Hole;
  G4SubtractionSolid* S800GateValve;

  // logical volumes
  G4LogicalVolume* S800_log;
  G4LogicalVolume* S800GateValve_log;
 
  // physical volumes
  G4VPhysicalVolume* S800_phys;
  G4VPhysicalVolume* S800GateValve_phys;

};

#endif
