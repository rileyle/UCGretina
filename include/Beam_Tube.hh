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

class Beam_Tube 
{
  public:

  G4LogicalVolume *expHall_log;

  Beam_Tube(G4LogicalVolume*,Materials*);
  ~Beam_Tube();
  
  G4VPhysicalVolume *Construct();
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

  //materials
  Materials* materials;
  G4Material* BeamTubeMaterial;

  //default position
  G4RotationMatrix NoRot;
  G4ThreeVector *Pos0;

  //the tube
  G4Tubs* BeamTube;

  //logical volume
  G4LogicalVolume* BeamTube_log;
 
  //physical volume
  G4VPhysicalVolume* BeamTube_phys;

};

#endif
#endif
