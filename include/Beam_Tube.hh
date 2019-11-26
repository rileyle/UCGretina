#ifndef LHTARGET
#ifndef Beam_Tube_H
#define Beam_Tube_H 1


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

class Beam_Tube 
{
  public:

  G4LogicalVolume *expHall_log;

  Beam_Tube(Materials*);
  ~Beam_Tube();
  
  G4VPhysicalVolume *Construct(G4LogicalVolume*);
  void setRmin(G4double);
  void setRmax(G4double);
  void setLength(G4double);
  void setMaterial(G4String);
  void Report();
  
    private:
  // dimensions
  G4double BTrmin;
  G4double BTrmax;
  G4double BTDz;
  G4double BTSPhi;
  G4double BTDPhi; 
  G4double BTOffset;
  G4double BTFlangeRmin;
  G4double BTFlangeRmax;
  G4double BTFlangeDz;
  G4double BTFlangeOffset;

  //materials
  Materials* materials;
  G4Material* BeamTubeMaterial;

  //default position
  G4RotationMatrix NoRot;
  G4ThreeVector *Pos0;
  G4ThreeVector *BTPos;
  G4ThreeVector *BTFlangePos;

  //the tube
  G4Tubs* BeamTube;
  G4Tubs* BeamTubeFlange;

  //logical volume
  G4LogicalVolume* BeamTube_log;
  G4LogicalVolume* BeamTubeFlange_log;
 
  //physical volume
  G4VPhysicalVolume* BeamTube_phys;
  G4VPhysicalVolume* BeamTubeFlange_phys;

};

#endif
#endif
