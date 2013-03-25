#ifndef Target_LH_H
#define Target_LH_H 1

#include "G4Material.hh"
#include "G4Element.hh"
#include "Materials.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Torus.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4AssemblyVolume.hh"

class Target 
{
  public:

  G4LogicalVolume *expHall_log;

  Target(G4LogicalVolume*,Materials*);
  ~Target();
  
  G4VPhysicalVolume *Construct();
  void setX(G4double);
  void setY(G4double);
  void setZ(G4double);
  void setMaterial(G4String);
  void setNStep(G4int);
  void setSourceFrame(G4String);
  void setSled();
  void Report();
  G4LogicalVolume* GetTargetLog(){return Target_log;}
  G4UnionSolid* GetTarget(){return aTarget;}  //was G4Tubs*
  G4VPhysicalVolume* GetTargetPlacement(){return Target_phys;}
  void setTargetReactionDepth(G4double);
  void ScaleDensity(G4double);
  void SetPositionZ(G4double);
 /* void setX(G4double);
  void setY(G4double);
  void setZ(G4double);*/
  G4double GetTargetThickness(){return Target_thickness;}
  G4ThreeVector* GetPos(){return Pos;}

private:
  // LH target beam tube dimensions
  G4double BTrmin;
  G4double BTrmax;
  G4double BTDz;
  G4double BTSPhi;
  G4double BTDPhi; 

  G4double LHTrmin;
  G4double LHTrmax;
  G4double LHTDz;
  G4double LHTSPhi;
  G4double LHTDPhi;

  G4double SubCylinderrmin;
  G4double SubCylinderrmax;
  G4double SubCylinderDz;
  G4double SubCylinderSPhi;
  G4double SubCylinderDPhi;

  // LH target dimensions
  G4double Target_rmax;
  G4double Target_rmin;
  G4double Target_z;
  G4double Target_thickness;
  G4double Target_SPhi;
  G4double Target_DPhi;
  G4double CnctrMaterial_gmol;
  G4double CnctrMaterial_density;
  G4double TargetCellMaterial_density;

  G4double BulgeLx_r;
  G4double BulgeLy_r;
  G4double BulgeLz_r;
  G4double BulgeL_cut;

  G4double BulgeRx_r;
  G4double BulgeRy_r;
  G4double BulgeRz_r;
  G4double BulgeR_cut;
  
  G4double TargetCell_rmax;
  G4double TargetCell_rmin;
  G4double TargetCell_z;
  G4double TargetCell_SPhi;
  G4double TargetCell_DPhi;
  
  G4double Cnctr1rmin;
  G4double Cnctr1rmax;
  G4double Cnctr1SPhi;
  G4double Cnctr1DPhi;
  G4double Cnctr1z;
  G4double Cnctr1y;

  G4double Cnctr2rmin;
  G4double Cnctr2rmax;
  G4double Cnctr2SPhi;
  G4double Cnctr2DPhi;
  G4double Cnctr2z;
  G4double Cnctr2y;

  G4double Cnctr3rmin;
  G4double Cnctr3rmax;
  G4double Cnctr3SPhi;
  G4double Cnctr3DPhi;
  G4double Cnctr3z;
  G4double Cnctr3y;

  G4double SteelTubermin;
  G4double SteelTubermax;
  G4double SteelTubeSPhi;
  G4double SteelTubeDPhi;
  G4double SteelTubez;
  G4double SteelTubey; 

  G4double Cnctr4stemrmin;
  G4double Cnctr4stemrmax;
  G4double Cnctr4stemSPhi;
  G4double Cnctr4stemDPhi;
  G4double Cnctr4stemz;
  G4double Cnctr4stemy;

  G4double Cnctr4rmin;
  G4double Cnctr4rmax;
  G4double Cnctr4SPhi;
  G4double Cnctr4DPhi;
  G4double Cnctr4z;
  G4double Cnctr4y;

  G4double CopperFlangermin;
  G4double CopperFlangermax;
  G4double CopperFlangeSPhi;
  G4double CopperFlangeDPhi;
  G4double CopperFlangez;
  G4double CopperFlangey;

  G4double SteelTorusrmin;
  G4double SteelTorusrmax;
  G4double SteelTorusrtor;
  G4double SteelTorusSPhi;
  G4double SteelTorusDPhi;
  G4double SteelTorusy;

  G4double Heaterrmin;
  G4double Heaterrmax;
  G4double HeaterSPhi;
  G4double HeaterDPhi;
  G4double Heaterz;
  G4double Heatery; 

  G4double secstageflangermin;
  G4double secstageflangermax;
  G4double secstageflangeSPhi;
  G4double secstageflangeDPhi;
  G4double secstageflangez;
  G4double secstageflangey;

  G4double secstage1rmin;
  G4double secstage1rmax;
  G4double secstage1SPhi;
  G4double secstage1DPhi;
  G4double secstage1z;
  G4double secstage1y;

  G4double secstage2rmin;
  G4double secstage2rmax;
  G4double secstage2SPhi;
  G4double secstage2DPhi;
  G4double secstage2z;
  G4double secstage2y;

  G4double fststageflangermin;
  G4double fststageflangermax;
  G4double fststageflangeSPhi;
  G4double fststageflangeDPhi;
  G4double fststageflangez;
  G4double fststageflangey;

  G4double fststage1rmin;
  G4double fststage1rmax;
  G4double fststage1SPhi;
  G4double fststage1DPhi;
  G4double fststage1z;
  G4double fststage1y;

  G4double fststageringrmin;
  G4double fststageringrmax;
  G4double fststageringSPhi;
  G4double fststageringDPhi;
  G4double fststageringz;
  G4double fststageringy;

  G4double fststage2rmin;
  G4double fststage2rmax;
  G4double fststage2SPhi;
  G4double fststage2DPhi;
  G4double fststage2z;
  G4double fststage2y;

  G4double coldhead1rmin;
  G4double coldhead1rmax;
  G4double coldhead1SPhi;
  G4double coldhead1DPhi;
  G4double coldhead1z;
  G4double coldhead1y;

  G4double coldhead2rmin;
  G4double coldhead2rmax;
  G4double coldhead2SPhi;
  G4double coldhead2DPhi;
  G4double coldhead2z;
  G4double coldhead2y;

  G4double coldhead3_x;
  G4double coldhead3_y;
  G4double coldhead3_z;
  G4double coldhead3z;
  G4double coldhead3y;

  G4double coldhead4rmin;
  G4double coldhead4rmax;
  G4double coldhead4SPhi;
  G4double coldhead4DPhi;
  G4double coldhead4z;
  G4double coldhead4y;

  G4double coldhead5rmin;
  G4double coldhead5rmax;
  G4double coldhead5SPhi;
  G4double coldhead5DPhi;
  G4double coldhead5z;
  G4double coldhead5y;

  G4double coldheadcylrmin;
  G4double coldheadcylrmax;
  G4double coldheadcylSPhi;
  G4double coldheadcylDPhi;
  G4double coldheadcylz;
  G4double coldheadcyly;
  G4double coldheadcylzpos;

  G4double alshieldrmin;
  G4double alshieldrmax;
  G4double alshieldSPhi;
  G4double alshieldDPhi;
  G4double alshieldz;
  G4double alshieldy;
  G4double cell_z;

  G4double bellowsdisc1rmin;
  G4double bellowsdisc1rmax;
  G4double bellowsdisc1SPhi;
  G4double bellowsdisc1DPhi;
  G4double bellowsdisc1z;
  G4double bellowsdisc1y;

  G4double bellowstorus1rmin;
  G4double bellowstorus1rmax;
  G4double bellowstorus1rtor;
  G4double bellowstorus1SPhi;
  G4double bellowstorus1DPhi;
  G4double bellowstorus1y;

  G4double bellowstorus2rmin;
  G4double bellowstorus2rmax;
  G4double bellowstorus2rtor;
  G4double bellowstorus2SPhi;
  G4double bellowstorus2DPhi;
  G4double bellowstorus2y;

  G4double bellowstorus3rmin;
  G4double bellowstorus3rmax;
  G4double bellowstorus3rtor;
  G4double bellowstorus3SPhi;
  G4double bellowstorus3DPhi;
  G4double bellowstorus3y;

  G4double bellowstorus4rmin;
  G4double bellowstorus4rmax;
  G4double bellowstorus4rtor;
  G4double bellowstorus4SPhi;
  G4double bellowstorus4DPhi;
  G4double bellowstorus4y;

  G4double bellowstorus5rmin;
  G4double bellowstorus5rmax;
  G4double bellowstorus5rtor;
  G4double bellowstorus5SPhi;
  G4double bellowstorus5DPhi;
  G4double bellowstorus5y;

  G4double bellowstorus6rmin;
  G4double bellowstorus6rmax;
  G4double bellowstorus6rtor;
  G4double bellowstorus6SPhi;
  G4double bellowstorus6DPhi;
  G4double bellowstorus6y;

  G4double bellowstorus7rmin;
  G4double bellowstorus7rmax;
  G4double bellowstorus7rtor;
  G4double bellowstorus7SPhi;
  G4double bellowstorus7DPhi;
  G4double bellowstorus7y;

  G4double bellowstorus8rmin;
  G4double bellowstorus8rmax;
  G4double bellowstorus8rtor;
  G4double bellowstorus8SPhi;
  G4double bellowstorus8DPhi;
  G4double bellowstorus8y;

  G4double bellowstorus9rmin;
  G4double bellowstorus9rmax;
  G4double bellowstorus9rtor;
  G4double bellowstorus9SPhi;
  G4double bellowstorus9DPhi;
  G4double bellowstorus9y;

  G4double bellowstorus10rmin;
  G4double bellowstorus10rmax;
  G4double bellowstorus10rtor;
  G4double bellowstorus10SPhi;
  G4double bellowstorus10DPhi;
  G4double bellowstorus10y;

  G4double bellowsdisc2rmin;
  G4double bellowsdisc2rmax;
  G4double bellowsdisc2z;
  G4double bellowsdisc2SPhi;
  G4double bellowsdisc2DPhi;
  G4double bellowsdisc2y;

  G4double sixptcnctrdisc1rmin;
  G4double sixptcnctrdisc1rmax;
  G4double sixptcnctrdisc1z;
  G4double sixptcnctrdisc1SPhi;
  G4double sixptcnctrdisc1DPhi;
  G4double sixptcnctrdisc1y;
 
  G4double sixptcnctrrmin;
  G4double sixptcnctrrmax;
  G4double sixptcnctrz;
  G4double sixptcnctrSPhi;
  G4double sixptcnctrDPhi;
  G4double sixptcnctry;

  G4double subcross1rmin;
  G4double subcross1rmax;
  G4double subcross1z;
  G4double subcross1SPhi;
  G4double subcross1DPhi;
  G4double subcross1y;

  G4double subcross2rmin;
  G4double subcross2rmax;
  G4double subcross2z;
  G4double subcross2SPhi;
  G4double subcross2DPhi;
  G4double subcross2y;

  G4double cross1rmin;
  G4double cross1rmax;
  G4double cross1z;
  G4double cross1SPhi;
  G4double cross1DPhi;
  G4double cross1y;
  
  G4double cross2rmin;
  G4double cross2rmax;
  G4double cross2z;
  G4double cross2SPhi;
  G4double cross2DPhi;
  G4double cross2y;

  G4double sixptcnctrdisc2rmin;
  G4double sixptcnctrdisc2rmax;
  G4double sixptcnctrdisc2z;
  G4double sixptcnctrdisc2SPhi;
  G4double sixptcnctrdisc2DPhi;
  G4double sixptcnctrdisc2y;

  G4double Steel_density;
  G4double lH2_density;
  G4double H_density;


  //source frame dimensions
  G4double frameThickness;
  G4double frameInnerRadius;
  G4double frameOuterRadius;
  G4double tapeThickness;
  G4double tape_r;
  
  //elements
  G4Element* elFe;
  G4Element* elC;
  G4Element* elCr;
  G4Element* elNi;
  G4Element* elMn;
  G4Element* elSi;
  G4Element* elP;
  G4Element* elS;
  G4Element* elH;
  
  //materials
  Materials* materials;
  G4Material* Final_BTMaterial;
  G4Material* TargetMaterial;
  G4Material* TargetCellMaterial;
  G4Material* Cnctr1Material;
  G4Material* Cnctr2Material;
  G4Material* Cnctr3Material;
  G4Material* Cnctr4stemMaterial;
  G4Material* Cnctr4Material;
  G4Material* CopperFlangeMaterial;
  G4Material* SteelTubeMaterial;
  G4Material* SteelTorusMaterial;
  G4Material* HeaterMaterial;
  G4Material* secstageflangeMaterial;
  G4Material* secstage1Material;
  G4Material* secstage2Material;
  G4Material* fststageflangeMaterial;
  G4Material* fststage1Material;
  G4Material* fststageringMaterial;
  G4Material* fststage2Material;
  G4Material* coldhead1Material;
  G4Material* coldhead2Material;
  G4Material* coldhead3Material;
  G4Material* coldhead4Material;
  G4Material* coldhead5Material;
  G4Material* alshieldMaterial;
  G4Material* frameMaterial;
  G4Material* tapeMaterial;
  G4Material* bellowsdisc1Material;
  G4Material* bellowstorus1Material;
  G4Material* bellowstorus2Material;
  G4Material* bellowstorus3Material;
  G4Material* bellowstorus4Material;
  G4Material* bellowstorus5Material;
  G4Material* bellowstorus6Material;
  G4Material* bellowstorus7Material;
  G4Material* bellowstorus8Material;
  G4Material* bellowstorus9Material;
  G4Material* bellowstorus10Material;
  G4Material* bellowsdisc2Material;
  G4Material* sixptcnctrdisc1Material;
  G4Material* sixptcnctrMaterial;
  G4Material* sixptcnctrdisc2Material;
  G4Material* cross1Material;
  G4Material* cross2Material;
  G4Material* coldheadcylMaterial;
  G4Material* Steel;
  G4Material* lH2;

  //default position
  G4RotationMatrix NoRot;
  G4ThreeVector *Pos;
  G4RotationMatrix LHTubeRot;
  G4RotationMatrix SubCylinderRot;
  G4RotationMatrix LHTargetRot;
  G4RotationMatrix Cnctr1Rot;
  G4RotationMatrix Cnctr4Rot;
  G4RotationMatrix subcross2Rot;
  G4ThreeVector *Cnctr1Pos;
  G4ThreeVector *Cnctr2Pos;
  G4ThreeVector *Cnctr3Pos;
  G4ThreeVector *Cnctr4stemPos;
  G4ThreeVector *Cnctr4Pos;
  G4ThreeVector *CopperFlangePos;
  G4ThreeVector *SteelTorusPos;
  G4ThreeVector *SteelTubePos;
  G4ThreeVector *HeaterPos;
  G4ThreeVector *secstageflangePos;
  G4ThreeVector *secstage1Pos;
  G4ThreeVector *secstage2Pos;
  G4ThreeVector *fststageflangePos;
  G4ThreeVector *fststage1Pos;
  G4ThreeVector *fststageringPos;
  G4ThreeVector *fststage2Pos;
  G4ThreeVector *coldhead1Pos;
  G4ThreeVector *coldhead2Pos;
  G4ThreeVector *coldhead3Pos;
  G4ThreeVector *coldhead4Pos;
  G4ThreeVector *coldhead5Pos;
  G4ThreeVector *alshieldPos;
  G4ThreeVector *BulgeLPos;
  G4ThreeVector *BulgeRPos;
  G4ThreeVector *bellowsdisc1Pos;
  G4ThreeVector *bellowstorus1Pos;
  G4ThreeVector *bellowstorus2Pos;
  G4ThreeVector *bellowstorus3Pos;
  G4ThreeVector *bellowstorus4Pos;
  G4ThreeVector *bellowstorus5Pos;
  G4ThreeVector *bellowstorus6Pos;
  G4ThreeVector *bellowstorus7Pos;
  G4ThreeVector *bellowstorus8Pos;
  G4ThreeVector *bellowstorus9Pos;
  G4ThreeVector *bellowstorus10Pos;
  G4ThreeVector *bellowsdisc2Pos;
  G4ThreeVector *sixptcnctrdisc1Pos;
  G4ThreeVector *sixptcnctrPos;
  G4ThreeVector *cross1Pos;
  G4ThreeVector *cross2Pos;
  G4ThreeVector *sixptcnctrdisc2Pos;
  G4ThreeVector *coldheadcylPos;
 
  //the shapes
  G4Tubs* BeamTube;
  G4Tubs* LHTargetTube; 
  G4Tubs* SubCylinder;

  G4Tubs* Target_;
  G4Tubs* aTargetCell;
  G4Tubs* Cnctr1;
  G4Tubs* Cnctr2;
  G4Tubs* Cnctr3;
  G4Tubs* Cnctr4stem;
  G4Tubs* Cnctr4;
  G4Tubs* CopperFlange;
  G4Tubs* SteelTube;
  G4Torus* SteelTorus;
  G4Tubs* Heater;
  G4Tubs* secstageflange;
  G4Tubs* secstage1;
  G4Tubs* secstage2;
  G4Tubs* fststageflange;
  G4Tubs* fststage1;
  G4Tubs* fststage2;
  G4Tubs* fststagering;
  G4Tubs* coldhead1;
  G4Tubs* coldhead2;
  G4Box* coldhead3;
  G4Tubs* coldhead4;
  G4Tubs* coldhead5;
  G4Tubs* coldheadcyl;
  G4Tubs* alshield;
  G4Tubs* cellsub;
  G4Ellipsoid* BulgeL;
  G4Ellipsoid* BulgeR;
  G4Tubs* euFrame;
  G4Tubs* euTape;
  G4Tubs* csFrame;
  G4Tubs* csRing;
  G4Tubs* csTape;
  G4Tubs* bellowsdisc1;
  G4Torus* bellowstorus1;
  G4Torus* bellowstorus2;
  G4Torus* bellowstorus3;
  G4Torus* bellowstorus4;
  G4Torus* bellowstorus5;
  G4Torus* bellowstorus6;
  G4Torus* bellowstorus7;
  G4Torus* bellowstorus8;
  G4Torus* bellowstorus9;
  G4Torus* bellowstorus10;
  G4Tubs* bellowsdisc2;
  G4Tubs* sixptcnctrdisc1;
  G4Tubs* sixptcnctr1;
  G4Tubs* subcross1;
  G4Tubs* subcross2;
  G4Tubs* cross1;
  G4Tubs* cross2;
  G4Tubs* sixptcnctrdisc2;

  //Assembly Volume
  G4AssemblyVolume* LHTarget;

  //Subtraction Solid
  G4SubtractionSolid* BeamTubeminusSubCylinder;
  G4SubtractionSolid* LHTargetTubeminusSubCylinder;
  G4SubtractionSolid* CnctrminusCell;
  G4SubtractionSolid* sixptcnctr2;
  G4SubtractionSolid* sixptcnctr3;
  G4SubtractionSolid* cross1_;
  G4SubtractionSolid* alshield_sub;

  //Union Solid
  G4UnionSolid* Final_BT;
  G4UnionSolid* Target_2;
  G4UnionSolid* aTarget;

  //logical volume
  G4LogicalVolume* Final_BT_log;
  G4LogicalVolume* Final_LHT_log;
  G4LogicalVolume* Target_log;
  G4LogicalVolume* TargetCell_log;
  G4LogicalVolume* Cnctr1_log;
  G4LogicalVolume* Cnctr2_log;
  G4LogicalVolume* Cnctr3_log;
  G4LogicalVolume* Cnctr4stem_log;
  G4LogicalVolume* Cnctr4_log;
  G4LogicalVolume* CopperFlange_log;
  G4LogicalVolume* SteelTorus_log;
  G4LogicalVolume* SteelTube_log;
  G4LogicalVolume* Heater_log;
  G4LogicalVolume* secstageflange_log;
  G4LogicalVolume* secstage1_log;
  G4LogicalVolume* secstage2_log;
  G4LogicalVolume* fststageflange_log;
  G4LogicalVolume* fststage1_log;
  G4LogicalVolume* fststagering_log;
  G4LogicalVolume* fststage2_log;
  G4LogicalVolume* coldhead1_log;
  G4LogicalVolume* coldhead2_log;
  G4LogicalVolume* coldhead3_log;
  G4LogicalVolume* coldhead4_log;
  G4LogicalVolume* coldhead5_log;
  G4LogicalVolume* coldheadcyl_log;
  G4LogicalVolume* alshield_log;
  G4LogicalVolume* euFrame_log;
  G4LogicalVolume* euTape_log;
  G4LogicalVolume* csFrame_log;
  G4LogicalVolume* csRing_log;
  G4LogicalVolume* csTape_log;
  G4LogicalVolume* bellowsdisc1_log;
  G4LogicalVolume* bellowstorus1_log;
  G4LogicalVolume* bellowstorus2_log;
  G4LogicalVolume* bellowstorus3_log;
  G4LogicalVolume* bellowstorus4_log;
  G4LogicalVolume* bellowstorus5_log;
  G4LogicalVolume* bellowstorus6_log;
  G4LogicalVolume* bellowstorus7_log;
  G4LogicalVolume* bellowstorus8_log;
  G4LogicalVolume* bellowstorus9_log;
  G4LogicalVolume* bellowstorus10_log;
  G4LogicalVolume* bellowsdisc2_log;
  G4LogicalVolume* sixptcnctrdisc1_log;
  G4LogicalVolume* sixptcnctr_log;
  G4LogicalVolume* cross1_log;
  G4LogicalVolume* cross2_log;
  G4LogicalVolume* sixptcnctrdisc2_log;

  
  //physical volume
  G4VPhysicalVolume* Target_phys;
  G4VPhysicalVolume* TargetCell_phys;
  G4VPhysicalVolume* Cnctr1_phys;
  G4VPhysicalVolume* Cnctr2_phys;
  G4VPhysicalVolume* Cnctr3_phys;
  G4VPhysicalVolume* Cnctr4stem_phys;
  G4VPhysicalVolume* SteelTorus_phys;
  G4VPhysicalVolume* CopperFlange_phys;
  G4VPhysicalVolume* Cnctr4_phys;
  G4VPhysicalVolume* SteelTube_phys;
  G4VPhysicalVolume* Heater_phys;
  G4VPhysicalVolume* secstageflange_phys;
  G4VPhysicalVolume* secstage1_phys;
  G4VPhysicalVolume* secstage2_phys;
  G4VPhysicalVolume* fststageflange_phys;
  G4VPhysicalVolume* fststage1_phys;
  G4VPhysicalVolume* fststage2_phys;
  G4VPhysicalVolume* fststagering_phys;
  G4VPhysicalVolume* coldhead1_phys;
  G4VPhysicalVolume* coldhead2_phys;
  G4VPhysicalVolume* coldhead3_phys;
  G4VPhysicalVolume* coldhead4_phys;
  G4VPhysicalVolume* coldhead5_phys;
  G4VPhysicalVolume* coldheadcyl_phys;
  G4VPhysicalVolume* alshield_phys;
  G4VPhysicalVolume* euFrame_phys;
  G4VPhysicalVolume* euTape_phys;
  G4VPhysicalVolume* csFrame_phys;
  G4VPhysicalVolume* csRing_phys;
  G4VPhysicalVolume* csTape_phys;
  G4VPhysicalVolume* bellowsdisc1_phys;
  G4VPhysicalVolume* bellowstorus1_phys;
  G4VPhysicalVolume* bellowstorus2_phys;
  G4VPhysicalVolume* bellowstorus3_phys;
  G4VPhysicalVolume* bellowstorus4_phys;
  G4VPhysicalVolume* bellowstorus5_phys;
  G4VPhysicalVolume* bellowstorus6_phys;
  G4VPhysicalVolume* bellowstorus7_phys;
  G4VPhysicalVolume* bellowstorus8_phys;
  G4VPhysicalVolume* bellowstorus9_phys;
  G4VPhysicalVolume* bellowstorus10_phys;
  G4VPhysicalVolume* bellowsdisc2_phys;
  G4VPhysicalVolume* sixptcnctrdisc1_phys;
  G4VPhysicalVolume* sixptcnctr_phys;
  G4VPhysicalVolume* cross1_phys;
  G4VPhysicalVolume* cross2_phys;
  G4VPhysicalVolume* sixptcnctrdisc2_phys;
 
  G4UserLimits *target_limits;

  //Number of simulation steps
  G4int NStep;

  G4String sourceFrame;

  G4double targetAngle;

};

#endif

