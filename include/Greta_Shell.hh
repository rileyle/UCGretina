////////////////////////////////////////////////////////////////////////////////
// This class provides the GRETINA mounting shell used at the NSCL. 
// (Lew Riley lriley@ursinus.edu) 
////////////////////////////////////////////////////////////////////////////////

#ifndef Greta_Shell_h
#define Greta_Shell_h 1

#include "DetectorConstruction.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4Material.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "globals.hh"

using namespace std;

class Greta_Shell
{
public:

  Greta_Shell();
  ~Greta_Shell();

private:

  /////////////////////////////
  /// Material and its name
  ////////////////////////////

  G4String     matShellName;
  G4Material  *matShell;
  
  /////////////////////////////
  /// Geometry parameters
  ////////////////////////////
  G4double         Rmin;
  G4double         Rmax;
  G4double         smallPortRadius;
  G4double         modulePortRadius;
  G4double         northOffset;
  G4double         southOffset;
  G4ThreeVector    Pos;
  G4ThreeVector    Pos0;
  G4RotationMatrix Rot;
  G4RotationMatrix Rot0;
  G4ThreeVector    PosSP[10];
  G4double         ModuleEuler[30][3];
  G4int            SmallPortStatus[10];
  G4int            ModulePortStatus[30];

public:
    G4int  FindMaterials();
    void   setNorthOffset(G4double off){northOffset = off;}
    void   setSouthOffset(G4double off){southOffset = off;}
    void   Placement(G4String);
    void   HalfShell(G4String);
};

#endif
