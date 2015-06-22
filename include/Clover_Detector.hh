#ifndef Clover_Detector_H
#define Clover_Detector_H 1

#include "G4Material.hh"
#include "Materials.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"

#include "TrackerGammaSD.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4AssemblyVolume.hh"

#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "Randomize.hh"
#include "globals.hh"
#include <iostream>
#include <iomanip>

using namespace std;

class Clover_Detector 
{
public:

  G4LogicalVolume *expHall_log;
  Materials* materials;
  
  Clover_Detector(G4LogicalVolume*, Materials*);
  ~Clover_Detector();

  G4VPhysicalVolume *Construct();

  void setX(G4double x);
  void setY(G4double y);
  void setZ(G4double z);

  void MakeSensitive(TrackerGammaSD*);

  private:

  // Logical volumes

  G4LogicalVolume* detector_log;
  G4LogicalVolume* fill_log;
  G4LogicalVolume* wall_log;
  G4LogicalVolume* corner_log;
  G4LogicalVolume* cover_log;
  G4LogicalVolume* Cubox_log;
  

  // Physical volumes
 
  G4VPhysicalVolume* detector_phys;
  G4VPhysicalVolume* wall_phys;
  G4VPhysicalVolume* wall1_phys;
  G4VPhysicalVolume* wall2_phys;
  G4VPhysicalVolume* wall3_phys;
  G4VPhysicalVolume* Leaf0_phys;
  G4VPhysicalVolume* Leaf1_phys;
  G4VPhysicalVolume* Leaf2_phys;
  G4VPhysicalVolume* Leaf3_phys;
  G4VPhysicalVolume* corner_phys;
  G4VPhysicalVolume* corner1_phys;
  G4VPhysicalVolume* corner2_phys;
  G4VPhysicalVolume* corner3_phys;
  G4VPhysicalVolume* cover_phys;
  G4VPhysicalVolume* Cubox_phys;
 
  // Materials
  G4Material* HpGe;
  G4Material* Al;
  G4Material* Cu;

  // dimensions
  G4double Length;
  G4double Radius;
  G4double startAngle;
  G4double spanningAngle;
  G4double CCoffset;
  G4double CCradius;
  G4double wallZoffset; 
  G4double walloffset;
  G4double corneroffset;
  G4double cornerRadius;
  G4double coverwidth;
  G4double coverlength;
  G4double coverthickness;
  G4double covergap;
  G4double Cuboxlength;
  G4double boxlength;
  G4double torusradius;

  //here
  G4Transform3D Dassembly;
  G4ThreeVector assemblypos;
  G4RotationMatrix assemblyrot;
  G4ThreeVector assemblyshift;
  G4AssemblyVolume* assemblyclover;

  // position
  G4RotationMatrix DetRot;
  G4RotationMatrix wallrot;
  G4ThreeVector detShift;
  G4ThreeVector Leaf1Shift;
  G4ThreeVector Leaf2Shift;
  G4ThreeVector Leaf3Shift;
  G4ThreeVector DetPos;
  G4ThreeVector Leaf1Pos;
  G4ThreeVector Leaf2Pos;
  G4ThreeVector Leaf3Pos;
  G4double thetad;
  G4double phid;
  G4double LeafShift;
  G4ThreeVector Leaf0Shift;
  G4ThreeVector Leaf0Pos;
  G4ThreeVector cornershift;
  G4ThreeVector cornerpos;
  G4ThreeVector corner1shift;
  G4ThreeVector corner1pos;
  G4ThreeVector corner2shift;
  G4ThreeVector corner2pos;
  G4ThreeVector corner3shift;
  G4ThreeVector corner3pos;
  G4ThreeVector covershift;
  G4ThreeVector coverpos;
  G4ThreeVector Cuboxshift;
  G4ThreeVector Cuboxpos;

  G4Tubs* detector;
  G4Tubs* CCsub;
  G4Tubs* fill;
  G4Torus* torus;
  G4Tubs* cornerCut;
  G4Tubs* corner;
  G4Tubs* CuboxCut;
  G4Tubs* torustube;
  G4Box* boxout;
  G4Box* boxin;
  G4Box* box;
  G4Box* Cubox;
  G4Box* torusbox;
  G4SubtractionSolid* cover;
  G4SubtractionSolid* coversub;
  G4SubtractionSolid* CuboxCut1;
  G4SubtractionSolid* CuboxCut2;
  G4SubtractionSolid* CuboxCut3;
  G4SubtractionSolid* CuboxCut4;
  G4SubtractionSolid* torus2;
  G4SubtractionSolid* detector_cut;
  G4UnionSolid* coveru;
  G4UnionSolid* torus1;
  G4UnionSolid* bevel;
  G4Box* wall;
  G4SubtractionSolid* subtract;
  G4IntersectionSolid* intersect;

};

#endif

