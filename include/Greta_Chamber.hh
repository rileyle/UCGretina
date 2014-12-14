#ifndef LHTARGET
#ifndef Greta_Chamber_H
#define Greta_Chamber_H 1


#include "G4Material.hh"
#include "Materials.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
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

class Greta_Chamber 
{
public:

  G4LogicalVolume *expHall_log;

  Greta_Chamber(G4LogicalVolume*,Materials*);
  ~Greta_Chamber();
  
  G4VPhysicalVolume *Construct();
  void setRmin(G4double r){ Rmin = r; }
  void setRmax(G4double r){ Rmax = r; }
  void setCutaway(){ Cutaway = true; }
  void Report();
  
private:
  // dimensions
  G4double Rmin;
  G4double Rmax;
  G4double BTrmin;
  G4double BTrmax;
  G4double BTDz;

  //materials
  Materials* materials;
  G4Material* ChamberMaterial;

  //logical volume
  G4LogicalVolume* Chamber_log;
 
  //physical volume
  G4VPhysicalVolume* Chamber_phys;

  G4bool Cutaway;

};

#endif
#endif
