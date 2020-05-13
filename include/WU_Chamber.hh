#ifndef LHTARGET
#ifndef WU_Chamber_H
#define WU_Chamber_H 1


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

class WU_Chamber 
{
  public:

  G4LogicalVolume *expHall_log;

  WU_Chamber();
  ~WU_Chamber();
  
  G4VPhysicalVolume *Construct(G4LogicalVolume*);
  void Report();
  
    private:
  // dimensions
  G4double BTrmin;
  G4double BTrmax;
  G4double BTDz;
  G4double BTSPhi;
  G4double BTDPhi; 
  G4double pRmin;
  G4double pRmax;
  G4double pSPhi;
  G4double pDPhi;
  G4double pSTheta;
  G4double pDTheta;
  G4double exBox;
  G4double smallCap;
  G4double cRmin;
  G4double cRmax;
  G4double boxSide;
  G4double boxSideLarge;
  G4double cosine1;
  G4double sine1;
  G4double apertureA;
  G4double apertureB;
  G4double apertureC;
  G4double RingDistance;
  G4double LTRadius1;
  G4double LTRadius2;
  G4double LTRadius3;
  G4double LTLength1;
  G4double LTLength2;
  G4double LTDistance;
  
  //materials
  G4Material* BeamTubeMaterial;
  G4Material* Brass;
  G4Material* Glass;

  //default position
  G4RotationMatrix NoRot;
  G4RotationMatrix CylRot;
  G4ThreeVector *Pos0;

  //the tube
  G4Tubs* BeamTube;
  G4Tubs* CylinderFill;

  //the hemisphere
  G4Sphere* Hemisphere;
  G4Sphere* Hemisphere2;

  //the connection tube
  G4Polycone* Connector;
  G4Polycone* Ring;
  G4Polycone* LongTube;
  G4Polycone* PortPlug;
  G4Polycone* TargetPort;

  //the other ports
  G4Tubs* Window;

  //logical volume
  G4LogicalVolume* BeamTube_log;
  G4LogicalVolume* Cap_log;
  G4LogicalVolume* CylinderFill_log;
  G4LogicalVolume* CylBox1_log;
  G4LogicalVolume* CylBox2_log;
  G4LogicalVolume* CylBox3_log;
  G4LogicalVolume* CylBox4_log;
  G4LogicalVolume* CylBox5_log;
  G4LogicalVolume* CylBox6_log;
  G4LogicalVolume* BeamHole1_log;
  G4LogicalVolume* BeamHole2_log;
  G4LogicalVolume* BeamHole3_log;
  G4LogicalVolume* BeamHole4_log;
  G4LogicalVolume* InCyl_log;
  G4LogicalVolume* Connector_log;
  G4LogicalVolume* Ring_log;
  G4LogicalVolume* LongTube_log;
  G4LogicalVolume* Window_log;
  G4LogicalVolume* PortPlug_log;
  G4LogicalVolume* Port_log;
  G4LogicalVolume* TargetPort_log;

  //physical volume
  G4VPhysicalVolume* BeamTube_phys;
  G4VPhysicalVolume* Cap_phys1;
  G4VPhysicalVolume* Cap_phys2;
  G4VPhysicalVolume* CylinderFill_phys;
  G4VPhysicalVolume* CylBox_phys;
  G4VPhysicalVolume* Connector_phys;
  G4VPhysicalVolume* Connector2_phys;
  G4VPhysicalVolume* Ring_phys;
  G4VPhysicalVolume* Ring2_phys;
  G4VPhysicalVolume* LongTube_phys;
  G4VPhysicalVolume* LongTube2_phys;
  G4VPhysicalVolume* Window_phys;
  G4VPhysicalVolume* PortPlug_phys;
  G4VPhysicalVolume* PortPlug2_phys;
  G4VPhysicalVolume* TargetPort_phys;

};

#endif
#endif
