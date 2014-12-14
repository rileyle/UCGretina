#ifndef Background_Sphere_H
#define Background_Sphere_H 1


#include "G4Material.hh"
#include "Materials.hh"
#include "G4Sphere.hh"
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

class Background_Sphere 
{
  public:

  G4LogicalVolume *expHall_log;

  Background_Sphere(G4LogicalVolume*,Materials*);
  ~Background_Sphere();
  
  G4VPhysicalVolume *Construct();
  void setRmin(G4double);
  void setRmax(G4double);
  void setMaterial(G4String);
  void Report();
  
  G4double GetRmin(){return BSrmin;}
  G4double GetRmax(){return BSrmax;}

    private:
  // dimensions
  G4double BSrmin;
  G4double BSrmax;

  //materials
  Materials* materials;
  G4Material* BackgroundSphereMaterial;

  //default position
  G4RotationMatrix NoRot;
  G4ThreeVector *Pos0;

  //the sphere
  G4Sphere* BackgroundSphere;

  //logical volume
  G4LogicalVolume* BackgroundSphere_log;
 
  //physical volume
  G4VPhysicalVolume* BackgroundSphere_phys;

};

#endif
