#ifdef LHTARGET
#include "Target_LH.hh"

Target::Target(G4LogicalVolume* experimentalHall_log,Materials* mat)
{
  materials=mat;
  expHall_log=experimentalHall_log;
  
  targetAngle = 30.0*deg;

  // Beam tube

  BTrmin=7.539*cm;
  BTrmax=7.639*cm;
  BTDz=44.*cm; //LR (approx target to gate valve flange)
  BTSPhi=0.*deg;
  BTDPhi=360*deg; 

  SubCylinderrmin=0*cm;
  SubCylinderrmax=7.62*cm;
  SubCylinderDz=18.5*cm;
  SubCylinderSPhi=0*deg;
  SubCylinderDPhi=360*deg;

  LHTrmin=4.992*cm;
  LHTrmax=5.092*cm;
  LHTDz=18.461*cm; 
  LHTSPhi=0.*deg;
  LHTDPhi=360*deg;
 
  Final_BTMaterial = materials->FindMaterial("Al");

  // LH target

  TargetCell_rmin=2.*cm;
  TargetCell_z=2.381*cm;

  Target_rmin=0.*cm;
  Target_rmax=TargetCell_rmin;
  Target_SPhi=0.*deg;
  Target_DPhi=360.*deg;
  Target_z=TargetCell_z;
  Pos = new G4ThreeVector(0.,0.,0.);

  //  H_density = 1.01*g/mole;
  //  lH2_density = .0678*g/cm3;
  //  elH = new G4Element("Hydrogen","H",1,H_density);
  //  lH2 = new G4Material("Liquid Hydrogen", lH2_density,1);
  //  lH2->AddElement(elH,2);
  
  TargetCell_rmax=2.699*cm;
  TargetCell_SPhi=0.*deg;
  TargetCell_DPhi=360.*deg;
  TargetCellMaterial = materials->FindMaterial("Al");
  elFe = new G4Element("Iron","Fe",26,55.845*g/mole);
  elC = new G4Element("Carbon","C",6,12.0107*g/mole);
  elCr = new G4Element("Chromium","Cr",24,51.9961*g/mole);
  elNi = new G4Element("Nickel","Ni",28,58.6934*g/mole);
  elMn = new G4Element("Manganese","Mn",25,54.938045*g/mole);
  elSi = new G4Element("Silicon","Si",14,28.0855*g/mole);
  elP = new G4Element("Phosphorus","P",15,30.973762*g/mole);
  elS = new G4Element("Sulfur","S",16,32.065*g/mole);
  Steel_density=.079*kg/cm3;
  Steel = new G4Material("Steel",Steel_density,8);
  Steel->AddElement(elFe,68.6*perCent);
  Steel->AddElement(elC,.08*perCent);
  Steel->AddElement(elCr,18.75*perCent);
  Steel->AddElement(elNi,9.5*perCent);
  Steel->AddElement(elMn,2*perCent);
  Steel->AddElement(elSi,1*perCent);
  Steel->AddElement(elP,.045*perCent);
  Steel->AddElement(elS,.025*perCent);

  Cnctr1rmin=0.*cm;
  Cnctr1rmax=.49*cm;
  Cnctr1SPhi=0.*deg;
  Cnctr1DPhi=360.*deg;
  Cnctr1z=(1.984/4)*cm; 
  Cnctr1y=TargetCell_rmax + Cnctr1z;
  Cnctr1Pos = new G4ThreeVector(0.,Cnctr1y,0.);
  Cnctr1Material = TargetCellMaterial;
  
  Cnctr2rmin=0.*cm;
  Cnctr2rmax=1.43*cm;
  Cnctr2SPhi=0.*deg;
  Cnctr2DPhi=360.*deg;
  Cnctr2z=(.635/2)*cm;
  Cnctr2y=Cnctr1y + (Cnctr1z)+(Cnctr2z);
  Cnctr2Pos= new G4ThreeVector(0.,Cnctr2y,0.);
  Cnctr2Material = TargetCellMaterial;
  
  Cnctr3rmin=0.*cm;
  Cnctr3rmax=1.43*cm;
  Cnctr3SPhi=0.*deg;
  Cnctr3DPhi=360.*deg;
  Cnctr3z=(.635/2)*cm;
  Cnctr3y=Cnctr2y + (Cnctr2z)+(Cnctr3z);
  Cnctr3Pos = new G4ThreeVector(0.,Cnctr3y,0.);
  Cnctr3Material = materials->FindMaterial("Cu");

  Cnctr4stemrmin=0.*cm;
  Cnctr4stemrmax=(.714/3)*cm;
  Cnctr4stemSPhi=0.*deg;
  Cnctr4stemDPhi=360.*deg;
  Cnctr4stemz=.5*cm;
  Cnctr4stemy=Cnctr3y+(Cnctr3z)+(Cnctr4stemz);
  Cnctr4stemPos = new G4ThreeVector(0.,Cnctr4stemy,0.);
  Cnctr4stemMaterial = materials->FindMaterial("Cu");

  Cnctr4rmin=0.*cm;
  Cnctr4rmax=.714*cm;
  Cnctr4SPhi=0.*deg;
  Cnctr4DPhi=360.*deg;
  Cnctr4z=(13.97/2)*cm;
  Cnctr4y=Cnctr4stemy + (Cnctr4stemz) + (Cnctr4z);
  Cnctr4Pos = new G4ThreeVector(0.,Cnctr4y,0.);
  Cnctr4Material = materials->FindMaterial("Cu");

  CopperFlangermin=0.*cm;
  CopperFlangermax=(2.619/2)*cm;
  CopperFlangeSPhi=0.*deg;
  CopperFlangeDPhi=360.*deg;
  CopperFlangez=(.635/2)*cm;
  CopperFlangey=Cnctr4y+Cnctr4z+.25*cm;
  CopperFlangePos = new G4ThreeVector(0.,CopperFlangey,0.);
  CopperFlangeMaterial = materials->FindMaterial("Cu");

  SteelTorusrmin=.159*cm;
  SteelTorusrmax=.318*cm;
  SteelTorusrtor= (3.175/4)*cm;
  SteelTorusSPhi=0.*deg;
  SteelTorusDPhi=270.*deg;
  SteelTorusy= Cnctr4stemy;
  SteelTorusPos= new G4ThreeVector(0.,SteelTorusy,0.);
  SteelTorusMaterial= Steel;
  
  SteelTubermin=.159*cm;
  SteelTubermax=.318*cm;
  SteelTubeSPhi=0.*deg;
  SteelTubeDPhi=360.*deg;
  SteelTubez=(26.3525/2)*cm; 
  SteelTubey= SteelTubez+Cnctr3y+SteelTorusrmax+.5*cm;
  SteelTubePos = new G4ThreeVector(.794*cm,SteelTubey,0.); 
  SteelTubeMaterial = Steel;

  Heaterrmin=0.*cm;
  Heaterrmax=(2.619/2)*cm;
  HeaterSPhi=0.*deg;
  HeaterDPhi=360.*deg;
  Heaterz= (.635/2)*cm;
  Heatery= CopperFlangey + CopperFlangez + Heaterz + (1/6)*cm;
  HeaterPos = new G4ThreeVector(0.,Heatery,0.);
  HeaterMaterial = materials->FindMaterial("Cu");

  secstageflangermin=0.*cm;
  secstageflangermax=(2.619/2)*cm;
  secstageflangeSPhi=0.*deg;
  secstageflangeDPhi=360.*deg;
  secstageflangez= (.3175/2)*cm;
  secstageflangey = Heatery + Heaterz + secstageflangez + (1/6)*cm;
  secstageflangePos = new G4ThreeVector(0.,secstageflangey,0.);
  secstageflangeMaterial = materials->FindMaterial("Cu");

  secstage1rmin=0.*cm;
  secstage1rmax=(1.746/2)*cm;
  secstage1SPhi=0.*deg;
  secstage1DPhi=360.*deg;
  secstage1z=(2.699/2)*cm;
  secstage1y= secstageflangey + secstageflangez + secstage1z + (1/6)*cm;
  secstage1Pos = new G4ThreeVector(0.,secstage1y,0.);
  secstage1Material = materials->FindMaterial("Cu");

  secstage2rmin=0.*cm;
  secstage2rmax=(1.429/2)*cm;
  secstage2SPhi=0.*deg;
  secstage2DPhi=360.*deg;
  secstage2z=(5.239/2)*cm;
  secstage2y= secstage1y + secstage1z + secstage2z + (1/6)*cm;
  secstage2Pos = new G4ThreeVector(0.,secstage2y,0.);
  secstage2Material = Steel;

  fststageflangermin = 0.*cm;
  fststageflangermax = (2.54/2)*cm;
  fststageflangeSPhi = 0.*deg;
  fststageflangeDPhi = 360.*deg;
  fststageflangez=(.159/2)*cm;
  fststageflangey = secstage2y + secstage2z + fststageflangez +(1/6)*cm;
  fststageflangePos = new G4ThreeVector(0.,fststageflangey,0.);
  fststageflangeMaterial = materials->FindMaterial("Cu");

  fststage1rmin = 0.*cm;
  fststage1rmax = (2.858/2)*cm;
  fststage1SPhi = 0.*deg;
  fststage1DPhi = 360.*deg;
  fststage1z= (2.858/2)*cm;
  fststage1y = fststageflangey + fststageflangez + fststage1z + (1/6)*cm;
  fststage1Pos = new G4ThreeVector(0.,fststage1y,0.);
  fststage1Material = materials->FindMaterial("Cu");

  fststageringrmin= fststage1rmax;
  fststageringrmax=(3.81/2)*cm;
  fststageringSPhi=0.*deg;
  fststageringDPhi=360.*deg;
  fststageringz=(.318/2)*cm;
  fststageringy = fststage1y + (fststageringz/2);
  fststageringPos = new G4ThreeVector(0.,fststageringy,0.);
  fststageringMaterial = materials->FindMaterial("Cu");

  fststage2rmin = 0.*cm;
  fststage2rmax = (2.699/2)*cm;
  fststage2SPhi = 0.*deg;
  fststage2DPhi = 360.*deg;
  fststage2z = (4.286/2)*cm;
  fststage2y = fststage1y + fststage1z + fststage2z + (1/6)*cm;
  fststage2Pos = new G4ThreeVector(0.,fststage2y,0.);
  fststage2Material = secstage2Material;

//want total bellowsz to be equal to 16.5*cm, max radii of tori= (2/3)*(16.5/10)*cm
  bellowsdisc1rmin=0.*cm;
  bellowsdisc1rmax=11.38*cm;
  bellowsdisc1SPhi=0.*deg;
  bellowsdisc1DPhi=360.*deg;
  bellowsdisc1z = (.318/2)*cm;
  bellowsdisc1y = fststage2y + fststage2z + bellowsdisc1z + (1/6)*cm;
  bellowsdisc1Pos = new G4ThreeVector(0.,bellowsdisc1y,0.);
  bellowsdisc1Material = secstage2Material;
  
  bellowstorus1rmin = 0.*cm;
  bellowstorus1rmax = 1.1*cm;
  bellowstorus1rtor = 7.639*cm;
  bellowstorus1SPhi = 0.*deg;
  bellowstorus1DPhi = 360.*deg;
  bellowstorus1y = bellowsdisc1y + bellowsdisc1z + bellowstorus1rmax + (1/6)*cm;
  bellowstorus1Pos = new G4ThreeVector(0.,bellowstorus1y,0.);
  bellowstorus1Material = secstage2Material;

  bellowstorus2rmin = 0.*cm;
  bellowstorus2rmax = 1.1*cm;
  bellowstorus2rtor = 7.639*cm;
  bellowstorus2SPhi = 0.*deg;
  bellowstorus2DPhi = 360.*deg;
  bellowstorus2y = bellowstorus1y + bellowstorus1rmax + bellowstorus2rmax + (1/6)*cm;
  bellowstorus2Pos = new G4ThreeVector(0.,bellowstorus2y,0.);
  bellowstorus2Material = secstage2Material;

  bellowstorus3rmin = 0.*cm;
  bellowstorus3rmax = 1.1*cm;
  bellowstorus3rtor = 7.639*cm;
  bellowstorus3SPhi = 0.*deg;
  bellowstorus3DPhi = 360.*deg;
  bellowstorus3y = bellowstorus2y + bellowstorus2rmax + bellowstorus3rmax + (1/6)*cm;
  bellowstorus3Pos = new G4ThreeVector(0.,bellowstorus3y,0.);
  bellowstorus3Material = secstage2Material;

  bellowstorus4rmin = 0.*cm;
  bellowstorus4rmax = 1.1*cm;
  bellowstorus4rtor = 7.639*cm;
  bellowstorus4SPhi = 0.*deg;
  bellowstorus4DPhi = 360.*deg;
  bellowstorus4y = bellowstorus3y + bellowstorus3rmax + bellowstorus4rmax + (1/6)*cm;
  bellowstorus4Pos = new G4ThreeVector(0.,bellowstorus4y,0.);
  bellowstorus4Material = secstage2Material;

  bellowstorus5rmin = 0.*cm;
  bellowstorus5rmax = 1.1*cm;
  bellowstorus5rtor = 7.639*cm;
  bellowstorus5SPhi = 0.*deg;
  bellowstorus5DPhi = 360.*deg;
  bellowstorus5y = bellowstorus4y + bellowstorus4rmax + bellowstorus5rmax + (1/6)*cm;
  bellowstorus5Pos = new G4ThreeVector(0.,bellowstorus5y,0.);
  bellowstorus5Material = secstage2Material;

  bellowstorus6rmin = 0.*cm;
  bellowstorus6rmax = 1.1*cm;
  bellowstorus6rtor = 7.639*cm;
  bellowstorus6SPhi = 0.*deg;
  bellowstorus6DPhi = 360.*deg;
  bellowstorus6y = bellowstorus5y + bellowstorus5rmax + bellowstorus6rmax + (1/6)*cm;
  bellowstorus6Pos = new G4ThreeVector(0.,bellowstorus6y,0.);
  bellowstorus6Material = secstage2Material;

  bellowstorus7rmin = 0.*cm;
  bellowstorus7rmax = 1.1*cm;
  bellowstorus7rtor = 7.639*cm;
  bellowstorus7SPhi = 0.*deg;
  bellowstorus7DPhi = 360.*deg;
  bellowstorus7y = bellowstorus6y + bellowstorus6rmax + bellowstorus7rmax + (1/6)*cm;
  bellowstorus7Pos = new G4ThreeVector(0.,bellowstorus7y,0.);
  bellowstorus7Material = secstage2Material;

  bellowstorus8rmin = 0.*cm;
  bellowstorus8rmax = 1.1*cm;
  bellowstorus8rtor = 7.639*cm;
  bellowstorus8SPhi = 0.*deg;
  bellowstorus8DPhi = 360.*deg;
  bellowstorus8y = bellowstorus7y + bellowstorus7rmax + bellowstorus8rmax + (1/6)*cm;
  bellowstorus8Pos = new G4ThreeVector(0.,bellowstorus8y,0.);
  bellowstorus8Material = secstage2Material;

  bellowstorus9rmin = 0.*cm;
  bellowstorus9rmax = 1.1*cm;
  bellowstorus9rtor = 7.639*cm;
  bellowstorus9SPhi = 0.*deg;
  bellowstorus9DPhi = 360.*deg;
  bellowstorus9y = bellowstorus8y + bellowstorus8rmax + bellowstorus9rmax + (1/6)*cm;
  bellowstorus9Pos = new G4ThreeVector(0.,bellowstorus9y,0.);
  bellowstorus9Material = secstage2Material;

  bellowstorus10rmin = 0.*cm;
  bellowstorus10rmax = 1.1*cm;
  bellowstorus10rtor = 7.639*cm;
  bellowstorus10SPhi = 0.*deg;
  bellowstorus10DPhi = 360.*deg;
  bellowstorus10y = bellowstorus9y + bellowstorus9rmax + bellowstorus10rmax + (1/6)*cm;
  bellowstorus10Pos = new G4ThreeVector(0.,bellowstorus10y,0.);
  bellowstorus10Material = secstage2Material;

  bellowsdisc2rmin = 7.366*cm;
  bellowsdisc2rmax = 11.38*cm;
  bellowsdisc2SPhi = 0.*deg;
  bellowsdisc2DPhi = 360.*deg;
  bellowsdisc2z = (.381/2)*cm;
  bellowsdisc2y = bellowstorus10y + bellowstorus10rmax + bellowsdisc2z + (1/6)*cm;
  bellowsdisc2Pos = new G4ThreeVector (0.,bellowsdisc2y,0.);
  bellowsdisc2Material = secstage2Material;

  sixptcnctrdisc1rmin = 7.366*cm;
  sixptcnctrdisc1rmax = 9.231*cm;
  sixptcnctrdisc1SPhi = 0.*deg;
  sixptcnctrdisc1DPhi = 360.*deg;
  sixptcnctrdisc1z = (.381/2)*cm;
  sixptcnctrdisc1y = bellowsdisc2y + bellowsdisc2z + sixptcnctrdisc1z + (1/6)*cm;
  sixptcnctrdisc1Pos = new G4ThreeVector (0.,sixptcnctrdisc1y,0.);
  sixptcnctrdisc1Material = Steel;

  sixptcnctrrmin = 7.366*cm;
  sixptcnctrrmax = 9.231*cm;
  sixptcnctrSPhi = 0.*deg;
  sixptcnctrDPhi = 360.*deg;
  sixptcnctrz = 5.08*cm;
  sixptcnctry = sixptcnctrdisc1y + sixptcnctrdisc1z + sixptcnctrz + (1/6)*cm;
  sixptcnctrPos = new G4ThreeVector (0.,sixptcnctry,0.);
  sixptcnctrMaterial = Steel;

  sixptcnctrdisc2rmin = 7.366*cm;
  sixptcnctrdisc2rmax = 7.639*cm;
  sixptcnctrdisc2SPhi = 0.*deg;
  sixptcnctrdisc2DPhi = 360.*deg;
  sixptcnctrdisc2z = (.381/2)*cm;
  sixptcnctrdisc2y = sixptcnctry + sixptcnctrz + sixptcnctrdisc2z + (1/6)*cm;
  sixptcnctrdisc2Pos = new G4ThreeVector (0.,sixptcnctrdisc2y,0.);
  sixptcnctrdisc2Material = Steel;

  subcross1rmin = 0.*cm;
  subcross1rmax = 2.54*cm;
  subcross1SPhi = 0.*deg;
  subcross1DPhi = 360.*deg; 
  subcross1z = 12.065*cm;
  subcross1y = sixptcnctry;
  
  subcross2rmin = 0.*cm;
  subcross2rmax = 2.54*cm;
  subcross2SPhi = 0.*deg;
  subcross2DPhi = 360.*deg;
  subcross2z = 12.065*cm;
  subcross2y = subcross1y;

  cross1rmin = 0.*cm;
  cross1rmax = 2.54*cm;
  cross1SPhi = 0.*deg;
  cross1DPhi = 360.*deg; 
  cross1z = subcross1z;
  cross1y = subcross1y;
  cross1Pos = new G4ThreeVector(0.,cross1y,0.);
  cross1Material = Steel; 
  
  cross2rmin = 0.*cm;
  cross2rmax = 2.54*cm;
  cross2SPhi = 0.*deg;
  cross2DPhi = 360.*deg;
  cross2z = subcross2z;
  cross2y = subcross2y;
  cross2Pos = new G4ThreeVector(0.,cross2y,0.);
  cross2Material = Steel;

  coldhead1rmin = 0.*cm;
  coldhead1rmax = 8.89*cm;
  coldhead1SPhi = 0.*deg;
  coldhead1DPhi = 360.*deg;
  coldhead1z = (.318/2)*cm;
  coldhead1y = sixptcnctrdisc2y + sixptcnctrdisc2z + (1/6)*cm;
  coldhead1Pos = new G4ThreeVector(0.,coldhead1y,0.);
  coldhead1Material = secstage2Material;

  coldhead2rmin = 0.*cm;
  coldhead2rmax = 7.639*cm;
  coldhead2SPhi = 0.*deg;
  coldhead2DPhi = 360.*deg;
  coldhead2z = (.476/2)*cm;
  coldhead2y = coldhead1y + coldhead1z + coldhead2z + (1/6)*cm;
  coldhead2Pos = new G4ThreeVector(0.,coldhead2y,0.);
  coldhead2Material = secstage2Material; 

  coldhead3_x = 6.056*cm;
  coldhead3_y = coldhead3_x;
  coldhead3_z = coldhead3_x;
  //coldhead3z = (6.033/2)*cm;
  coldhead3y = coldhead2y + coldhead2z + coldhead3_z + (1/6)*cm;
  coldhead3Pos = new G4ThreeVector(0.,coldhead3y,0.);
  coldhead3Material = secstage2Material;

  coldhead4rmin = 0.*cm;
  coldhead4rmax = 4.604*cm;
  coldhead4SPhi = 0.*deg;
  coldhead4DPhi = 360.*deg;
  coldhead4z = (.318/2)*cm;
  coldhead4y = coldhead3y + coldhead3_z + coldhead4z + (1/6)*cm;
  coldhead4Pos = new G4ThreeVector(0.,coldhead4y,0.);
  coldhead4Material = secstage2Material;

  coldhead5rmin = 0.*cm;
  coldhead5rmax = 1.429*cm;
  coldhead5SPhi = 0.*deg;
  coldhead5DPhi = 360.*deg;
  coldhead5z = (.953/2)*cm;
  coldhead5y = coldhead4y + coldhead4z + coldhead5z + (1/6)*cm;
  coldhead5Pos = new G4ThreeVector(0.,coldhead5y,0.);
  coldhead5Material = secstage2Material;

  coldheadcylrmin = 0.*cm;
  coldheadcylrmax = 4.542*cm;
  coldheadcylSPhi = 0.*deg;
  coldheadcylDPhi = 360.*deg;
  coldheadcylz = 5.57*cm;
  coldheadcyly = coldhead3y;
  coldheadcylzpos = -10.598*cm;
  coldheadcylPos = new G4ThreeVector(0.,coldheadcyly,coldheadcylzpos);
  coldheadcylMaterial = secstage2Material;

  alshieldrmin = 3.71*cm;
  alshieldrmax = 3.81*cm;
  alshieldSPhi = 0.*deg;
  alshieldDPhi = 360.*deg;
  alshieldz = 25.383*cm; 
  alshieldy = (20.383)*cm; 
  alshieldPos = new G4ThreeVector(0.,alshieldy,0.);
  alshieldMaterial = materials->FindMaterial("Al");
  cell_z = 5.*cm;

  BulgeLx_r=1.9*cm;
  BulgeLy_r=1.9*cm;
  BulgeLz_r=.1*cm;
  BulgeL_cut=0.*cm;
  BulgeLPos = new G4ThreeVector(0.,0.,Target_z/2.);

  BulgeRx_r=1.9*cm;
  BulgeRy_r=1.9*cm;
  BulgeRz_r=.1*cm;
  BulgeR_cut=0.*cm;
  BulgeRPos = new G4ThreeVector(0.,0.,-Target_z/2.);

  Target_thickness = Target_z + BulgeRz_r + BulgeLz_r; //to add LH2 back, add BulgeRz_r and BulgeLz_r

  lH2 = materials->FindMaterial("G4_lH2");

  NStep=20;
  
  sourceFrame = "";

}

Target::~Target()
{ delete target_limits;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* Target::Construct()
{
  LHTarget = new G4AssemblyVolume();

  // Beam tube

  
  BeamTube = new G4Tubs("BeamTube",BTrmin,BTrmax,BTDz,BTSPhi,BTDPhi); 
  SubCylinder = new G4Tubs("SubCylinder",SubCylinderrmin,SubCylinderrmax,SubCylinderDz,SubCylinderSPhi,SubCylinderDPhi);
  
  SubCylinderRot = G4RotationMatrix::IDENTITY;
  SubCylinderRot.rotateX(90.0*deg);

  NoRot = G4RotationMatrix::IDENTITY;

  BeamTubeminusSubCylinder = new G4SubtractionSolid ("Final_BT",BeamTube,SubCylinder,G4Transform3D(NoRot,G4ThreeVector(0.,LHTDz,0.)));
  
  LHTargetTube = new G4Tubs("LHTargetTube",LHTrmin,LHTrmax,LHTDz,LHTSPhi,LHTDPhi);

  LHTargetTubeminusSubCylinder = new G4SubtractionSolid ("Final_LHT",LHTargetTube,SubCylinder,G4Transform3D(SubCylinderRot,G4ThreeVector(0.,0.,LHTDz)));
 
  Final_BT = new G4UnionSolid("Final_BT", BeamTubeminusSubCylinder,LHTargetTubeminusSubCylinder,G4Transform3D(SubCylinderRot,G4ThreeVector(0.,LHTDz,0.))); 

  Final_BT_log = new G4LogicalVolume(Final_BT,Final_BTMaterial,"FinalBT_log",0,0,0);

  LHTarget->AddPlacedVolume(Final_BT_log, *Pos, &NoRot);

  // LH target

  Target_ = new G4Tubs("target",Target_rmin,Target_rmax,Target_z/2.,Target_SPhi,Target_DPhi);
  BulgeR = new G4Ellipsoid("bulger",BulgeRx_r,BulgeRy_r,BulgeRz_r);
  BulgeL = new G4Ellipsoid("bulgel",BulgeLx_r,BulgeLy_r,BulgeLz_r);
  
  Target_2 = new G4UnionSolid("target_2",Target_,BulgeR,G4Transform3D(NoRot,*BulgeRPos));
  aTarget = new G4UnionSolid("atarget",Target_2,BulgeL,G4Transform3D(NoRot,*BulgeLPos)); 

  Target_log = new G4LogicalVolume(aTarget,lH2,"target_log",0,0,0);
  target_limits= new G4UserLimits();
  target_limits->SetMaxAllowedStep(Target_thickness/NStep);
  Target_log->SetUserLimits(target_limits);
  Target_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),Target_log,"target",expHall_log,false,0);
  //  LHTarget->AddPlacedVolume(Target_log,*Pos, &NoRot);

  aTargetCell = new G4Tubs("TargetCell",TargetCell_rmin,TargetCell_rmax,TargetCell_z/2.,TargetCell_SPhi,TargetCell_DPhi);
  TargetCell_log = new G4LogicalVolume(aTargetCell,TargetCellMaterial,"TargetCell_log",0,0,0);
  //  TargetCell_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),TargetCell_log,"TargetCell",expHall_log,false,0);
  LHTarget->AddPlacedVolume(TargetCell_log,*Pos, &NoRot);

  Cnctr1Rot = G4RotationMatrix::IDENTITY;
  Cnctr1Rot.rotateX(90.0*deg);

  Cnctr1 = new G4Tubs("alcnctr",Cnctr1rmin,Cnctr1rmax,Cnctr1z,Cnctr1SPhi,Cnctr1DPhi);
  CnctrminusCell = new G4SubtractionSolid("CnctrminusCell",Cnctr1,aTargetCell,G4Transform3D(NoRot,G4ThreeVector(0.,4.55*cm,0.)));
  Cnctr1_log = new G4LogicalVolume(CnctrminusCell,Cnctr1Material,"Cnctr1_log",0,0,0);
  //  Cnctr1_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*Cnctr1Pos),Cnctr1_log,"Cnctr1",expHall_log,false,0);
  LHTarget->AddPlacedVolume(Cnctr1_log, *Cnctr1Pos, &Cnctr1Rot);

  Cnctr2 = new G4Tubs("alcnctr2",Cnctr2rmin,Cnctr2rmax,Cnctr2z,Cnctr2SPhi,Cnctr2DPhi);
  Cnctr2_log = new G4LogicalVolume(Cnctr2,Cnctr2Material,"Cnctr2_log",0,0,0);
  //  Cnctr2_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*Cnctr2Pos),Cnctr2_log,"Cnctr2",expHall_log,false,0);
  LHTarget->AddPlacedVolume(Cnctr2_log, *Cnctr2Pos, &Cnctr1Rot);

  Cnctr3 = new G4Tubs("cucnctr3",Cnctr3rmin,Cnctr3rmax,Cnctr3z,Cnctr3SPhi,Cnctr3DPhi);
  Cnctr3_log = new G4LogicalVolume(Cnctr3,Cnctr3Material,"Cnctr3_log",0,0,0);
  //Cnctr3_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*Cnctr3Pos),Cnctr3_log,"Cnctr3",expHall_log,false,0);
  LHTarget->AddPlacedVolume(Cnctr3_log, *Cnctr3Pos, &Cnctr1Rot);

  Cnctr4stem = new G4Tubs("cnctr4stem",Cnctr4stemrmin,Cnctr4stemrmax,Cnctr4stemz,Cnctr4stemSPhi,Cnctr4stemDPhi);
  Cnctr4stem_log = new G4LogicalVolume(Cnctr4stem,Cnctr4stemMaterial,"Cnctr4stem_log",0,0,0);
  //Cnctr4stem_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*Cnctr4stemPos),Cnctr4stem_log,"Cnctr4stem",expHall_log,false,0);
  LHTarget->AddPlacedVolume(Cnctr4stem_log, *Cnctr4stemPos, &Cnctr1Rot);

  Cnctr4 = new G4Tubs("cucnctr4",Cnctr4rmin,Cnctr4rmax,Cnctr4z,Cnctr4SPhi,Cnctr4DPhi);
  Cnctr4_log = new G4LogicalVolume(Cnctr4,Cnctr4Material,"Cnctr4_log",0,0,0);
  //  Cnctr4_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*Cnctr4Pos),Cnctr4_log,"Cnctr4",expHall_log,false,0);
  LHTarget->AddPlacedVolume(Cnctr4_log, *Cnctr4Pos, &Cnctr1Rot);

  CopperFlange = new G4Tubs("cuflange",CopperFlangermin,CopperFlangermax,CopperFlangez,CopperFlangeSPhi,CopperFlangeDPhi);
  CopperFlange_log = new G4LogicalVolume(CopperFlange,CopperFlangeMaterial,"CopperFlange_log",0,0,0);
  //  CopperFlange_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*CopperFlangePos),CopperFlange_log,"CopperFlange",expHall_log,false,0);
  LHTarget->AddPlacedVolume(CopperFlange_log, *CopperFlangePos, &Cnctr1Rot);

  SteelTorus = new G4Torus("steeltorus",SteelTorusrmin,SteelTorusrmax,SteelTorusrtor,SteelTorusSPhi,SteelTorusDPhi);
  SteelTorus_log = new G4LogicalVolume(SteelTorus,SteelTorusMaterial,"SteelTorus_log",0,0,0);
  //  SteelTorus_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*SteelTorusPos),SteelTorus_log,"SteelTorus",expHall_log,false,0);
  LHTarget->AddPlacedVolume(SteelTorus_log, *SteelTorusPos, &Cnctr1Rot);

  SteelTube = new G4Tubs("steeltube",SteelTubermin,SteelTubermax,SteelTubez,SteelTubeSPhi,SteelTubeDPhi);
  SteelTube_log = new G4LogicalVolume(SteelTube,SteelTubeMaterial,"SteelTube_log",0,0,0);
  //  SteelTube_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*SteelTubePos),SteelTube_log,"SteelTube",expHall_log,false,0);
  LHTarget->AddPlacedVolume(SteelTube_log, *SteelTubePos, &Cnctr1Rot);

  Heater = new G4Tubs("heater",Heaterrmin,Heaterrmax,Heaterz,HeaterSPhi,HeaterDPhi);
  Heater_log = new G4LogicalVolume(Heater,HeaterMaterial,"Heater_log",0,0,0);
  //  Heater_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*HeaterPos),Heater_log,"Heater",expHall_log,false,0);
  LHTarget->AddPlacedVolume(Heater_log, *HeaterPos, &Cnctr1Rot);

  secstageflange = new G4Tubs("secstageflange", secstageflangermin, secstageflangermax,secstageflangez,secstageflangeSPhi,secstageflangeDPhi);
  secstageflange_log = new G4LogicalVolume(secstageflange,secstageflangeMaterial,"secstageflange_log",0,0,0);
  //  secstageflange_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*secstageflangePos),secstageflange_log,"secstageflange",expHall_log,false,0);
  LHTarget->AddPlacedVolume(secstageflange_log, *secstageflangePos, &Cnctr1Rot);

  secstage1 = new G4Tubs("secstagept1",secstage1rmin,secstage1rmax,secstage1z,secstage1SPhi,secstage1DPhi);
  secstage1_log = new G4LogicalVolume(secstage1,secstage1Material,"secstage1_log",0,0,0);
  //  secstage1_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*secstage1Pos),secstage1_log,"secstage1",expHall_log,false,0);
  LHTarget->AddPlacedVolume(secstage1_log, *secstage1Pos, &Cnctr1Rot);

  secstage2 = new G4Tubs("secstagept2",secstage2rmin,secstage2rmax,secstage2z,secstage2SPhi,secstage2DPhi);
  secstage2_log = new G4LogicalVolume(secstage2,secstage2Material,"secstage2_log",0,0,0);
  //  secstage2_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*secstage2Pos),secstage2_log,"secstage2",expHall_log,false,0);
  LHTarget->AddPlacedVolume(secstage2_log, *secstage2Pos, &Cnctr1Rot);

  fststageflange = new G4Tubs("fststageflange",fststageflangermin,fststageflangermax,fststageflangez,fststageflangeSPhi,fststageflangeDPhi);
  fststageflange_log = new G4LogicalVolume(fststageflange,fststageflangeMaterial,"fststageflange_log",0,0,0);
  //  fststageflange_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*fststageflangePos),fststageflange_log,"fststageflange",expHall_log,false,0);
  LHTarget->AddPlacedVolume(fststageflange_log, *fststageflangePos, &Cnctr1Rot);

  fststage1 = new G4Tubs("fststagept1",fststage1rmin,fststage1rmax,fststage1z,fststage1SPhi,fststage1DPhi);
  fststage1_log = new G4LogicalVolume(fststage1,fststage1Material,"fststage1_log",0,0,0);
  //  fststage1_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*fststage1Pos),fststage1_log,"fststage1",expHall_log,false,0);
  LHTarget->AddPlacedVolume(fststage1_log, *fststage1Pos, &Cnctr1Rot);

  fststagering = new G4Tubs("ring",fststageringrmin,fststageringrmax,fststageringz,fststageringSPhi,fststageringDPhi);
  fststagering_log = new G4LogicalVolume(fststagering,fststageringMaterial,"fststagering_log",0,0,0);
  //  fststagering_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*fststageringPos),fststagering_log,"fststagering",expHall_log,false,0);
  LHTarget->AddPlacedVolume(fststagering_log, *fststageringPos, &Cnctr1Rot);

  fststage2 = new G4Tubs("fststagept2",fststage2rmin,fststage2rmax,fststage2z,fststage2SPhi,fststage2DPhi);
  fststage2_log = new G4LogicalVolume(fststage2,fststage2Material,"fststage2_log",0,0,0);
  //  fststage2_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*fststage2Pos),fststage2_log,"fststage2",expHall_log,false,0);
  LHTarget->AddPlacedVolume(fststage2_log, *fststage2Pos, &Cnctr1Rot);

  coldhead1 = new G4Tubs("coldheadpt1",coldhead1rmin,coldhead1rmax,coldhead1z,coldhead1SPhi,coldhead1DPhi);
  coldhead1_log = new G4LogicalVolume(coldhead1,coldhead1Material,"coldhead1_log",0,0,0);
  //  coldhead1_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*coldhead1Pos),coldhead1_log,"coldhead1",expHall_log,false,0);
  LHTarget->AddPlacedVolume(coldhead1_log, *coldhead1Pos, &Cnctr1Rot);

  coldhead2 = new G4Tubs("coldheadpt2",coldhead2rmin,coldhead2rmax,coldhead2z,coldhead2SPhi,coldhead2DPhi);
  coldhead2_log = new G4LogicalVolume(coldhead2,coldhead2Material,"coldhead2_log",0,0,0);
  //  coldhead2_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*coldhead2Pos),coldhead2_log,"coldhead2",expHall_log,false,0);
  LHTarget->AddPlacedVolume(coldhead2_log, *coldhead2Pos, &Cnctr1Rot);

  coldhead3 = new G4Box("coldheadpt3",coldhead3_x,coldhead3_y,coldhead3_z);
  coldhead3_log = new G4LogicalVolume(coldhead3,coldhead3Material,"coldhead3_log",0,0,0);
  //  coldhead3_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*coldhead3Pos),coldhead3_log,"coldhead3",expHall_log,false,0);
  LHTarget->AddPlacedVolume(coldhead3_log, *coldhead3Pos, &Cnctr1Rot);

  coldhead4 = new G4Tubs("coldheadpt4",coldhead4rmin,coldhead4rmax,coldhead4z,coldhead4SPhi,coldhead4DPhi);
  coldhead4_log = new G4LogicalVolume(coldhead4,coldhead4Material,"coldhead4_log",0,0,0);
  //  coldhead4_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*coldhead4Pos),coldhead4_log,"coldhead4",expHall_log,false,0);
  LHTarget->AddPlacedVolume(coldhead4_log, *coldhead4Pos, &Cnctr1Rot);

  coldhead5 = new G4Tubs("coldheadpt5",coldhead5rmin,coldhead5rmax,coldhead5z,coldhead5SPhi,coldhead5DPhi);
  coldhead5_log = new G4LogicalVolume(coldhead5,coldhead5Material,"coldhead5_log",0,0,0);
  //  coldhead5_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*coldhead5Pos),coldhead5_log,"coldhead5",expHall_log,false,0);
  LHTarget->AddPlacedVolume(coldhead5_log, *coldhead5Pos, &Cnctr1Rot);

  coldheadcyl = new G4Tubs("coldheadcyl",coldheadcylrmin,coldheadcylrmax,coldheadcylz,coldheadcylSPhi,coldheadcylDPhi);
  coldheadcyl_log = new G4LogicalVolume(coldheadcyl,coldheadcylMaterial,"coldheadcyl_log",0,0,0);
  //  coldheadcyl_phys = new G4PVPlacement(G4Transform3D(NoRot,*coldheadcylPos),coldheadcyl_log,"coldheadcyl",expHall_log,false,0);
  LHTarget->AddPlacedVolume(coldheadcyl_log, *coldheadcylPos, &NoRot);

  alshield = new G4Tubs("alshield",alshieldrmin,alshieldrmax,alshieldz,alshieldSPhi,alshieldDPhi);
  cellsub = new G4Tubs("cellsub",Target_rmin,TargetCell_rmax,cell_z,TargetCell_SPhi,TargetCell_DPhi);
  alshield_sub = new G4SubtractionSolid("alshieldsub",alshield,cellsub,G4Transform3D(Cnctr1Rot, G4ThreeVector(0.,0.,((alshieldy/2) + 10*cm))));
  alshield_log = new G4LogicalVolume(alshield_sub,alshieldMaterial,"alshield_log",0,0,0);
  //  alshield_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*alshieldPos),alshield_log,"alshield_sub",expHall_log,false,0);
  LHTarget->AddPlacedVolume(alshield_log, *alshieldPos, &Cnctr1Rot);

  bellowsdisc1 = new G4Tubs("bellowstopdisc",bellowsdisc1rmin,bellowsdisc1rmax,bellowsdisc1z,bellowsdisc1SPhi,bellowsdisc1DPhi);
  bellowsdisc1_log = new G4LogicalVolume(bellowsdisc1,bellowsdisc1Material,"bellowsdisc1_log",0,0,0);
  //  bellowsdisc1_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowsdisc1Pos),bellowsdisc1_log,"bellowsdisc1",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowsdisc1_log, *bellowsdisc1Pos, &Cnctr1Rot);

  bellowstorus1 = new G4Torus("torus1",bellowstorus1rmin,bellowstorus1rmax,bellowstorus1rtor,bellowstorus1SPhi,bellowstorus1DPhi);
  bellowstorus1_log = new G4LogicalVolume(bellowstorus1,bellowstorus1Material,"bellowstorus1_log",0,0,0);
  //  bellowstorus1_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus1Pos),bellowstorus1_log,"bellowstorus1",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus1_log, *bellowstorus1Pos, &Cnctr1Rot);

  bellowstorus2 = new G4Torus("torus2",bellowstorus2rmin,bellowstorus2rmax,bellowstorus2rtor,bellowstorus2SPhi,bellowstorus2DPhi);
  bellowstorus2_log = new G4LogicalVolume(bellowstorus2,bellowstorus2Material,"bellowstorus2_log",0,0,0);
  //  bellowstorus2_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus2Pos),bellowstorus2_log,"bellowstorus2",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus2_log, *bellowstorus2Pos, &Cnctr1Rot);

  bellowstorus3 = new G4Torus("torus3",bellowstorus3rmin,bellowstorus3rmax,bellowstorus3rtor,bellowstorus3SPhi,bellowstorus3DPhi);
  bellowstorus3_log = new G4LogicalVolume(bellowstorus3,bellowstorus3Material,"bellowstorus3_log",0,0,0);
  //  bellowstorus3_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus3Pos),bellowstorus3_log,"bellowstorus3",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus3_log, *bellowstorus3Pos, &Cnctr1Rot);

  bellowstorus4 = new G4Torus("torus4",bellowstorus4rmin,bellowstorus4rmax,bellowstorus4rtor,bellowstorus4SPhi,bellowstorus4DPhi);
  bellowstorus4_log = new G4LogicalVolume(bellowstorus4,bellowstorus4Material,"bellowstorus4_log",0,0,0);
  //  bellowstorus4_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus4Pos),bellowstorus4_log,"bellowstorus4",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus4_log, *bellowstorus4Pos, &Cnctr1Rot);

  bellowstorus5 = new G4Torus("torus5",bellowstorus5rmin,bellowstorus5rmax,bellowstorus5rtor,bellowstorus5SPhi,bellowstorus5DPhi);
  bellowstorus5_log = new G4LogicalVolume(bellowstorus5,bellowstorus5Material,"bellowstorus5_log",0,0,0);
  //  bellowstorus5_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus5Pos),bellowstorus5_log,"bellowstorus5",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus5_log, *bellowstorus5Pos, &Cnctr1Rot);

  bellowstorus6 = new G4Torus("torus6",bellowstorus6rmin,bellowstorus6rmax,bellowstorus6rtor,bellowstorus6SPhi,bellowstorus6DPhi);
  bellowstorus6_log = new G4LogicalVolume(bellowstorus6,bellowstorus6Material,"bellowstorus6_log",0,0,0);
  //  bellowstorus6_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus6Pos),bellowstorus6_log,"bellowstorus6",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus6_log, *bellowstorus6Pos, &Cnctr1Rot);

  bellowstorus7 = new G4Torus("torus7",bellowstorus7rmin,bellowstorus7rmax,bellowstorus7rtor,bellowstorus7SPhi,bellowstorus7DPhi);
  bellowstorus7_log = new G4LogicalVolume(bellowstorus7,bellowstorus7Material,"bellowstorus7_log",0,0,0);
  //  bellowstorus7_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus7Pos),bellowstorus7_log,"bellowstorus7",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus7_log, *bellowstorus7Pos, &Cnctr1Rot);

  bellowstorus8 = new G4Torus("torus8",bellowstorus8rmin,bellowstorus8rmax,bellowstorus8rtor,bellowstorus8SPhi,bellowstorus8DPhi);
  bellowstorus8_log = new G4LogicalVolume(bellowstorus8,bellowstorus8Material,"bellowstorus8_log",0,0,0);
  //  bellowstorus8_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus8Pos),bellowstorus8_log,"bellowstorus8",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus8_log, *bellowstorus8Pos, &Cnctr1Rot);

  bellowstorus9 = new G4Torus("torus9",bellowstorus9rmin,bellowstorus9rmax,bellowstorus9rtor,bellowstorus9SPhi,bellowstorus9DPhi);
  bellowstorus9_log = new G4LogicalVolume(bellowstorus9,bellowstorus9Material,"bellowstorus9_log",0,0,0);
  //  bellowstorus9_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus9Pos),bellowstorus9_log,"bellowstorus9",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus9_log, *bellowstorus9Pos, &Cnctr1Rot);

  bellowstorus10 = new G4Torus("torus10",bellowstorus10rmin,bellowstorus10rmax,bellowstorus10rtor,bellowstorus10SPhi,bellowstorus10DPhi);
  bellowstorus10_log = new G4LogicalVolume(bellowstorus10,bellowstorus10Material,"bellowstorus10_log",0,0,0);
  //  bellowstorus10_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowstorus10Pos),bellowstorus10_log,"bellowstorus10",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowstorus10_log, *bellowstorus10Pos, &Cnctr1Rot);

  bellowsdisc2 = new G4Tubs("bellowsdisc2",bellowsdisc2rmin,bellowsdisc2rmax,bellowsdisc2z,bellowsdisc2SPhi,bellowsdisc2DPhi);
  bellowsdisc2_log = new G4LogicalVolume(bellowsdisc2,bellowsdisc2Material,"bellowsdisc2_log",0,0,0);
  //  bellowsdisc2_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*bellowsdisc2Pos),bellowsdisc2_log,"bellowsdisc2",expHall_log,false,0);
  LHTarget->AddPlacedVolume(bellowsdisc2_log, *bellowsdisc2Pos, &Cnctr1Rot);

  sixptcnctrdisc1 = new G4Tubs("sixptcnctrdisc1",sixptcnctrdisc1rmin,sixptcnctrdisc1rmax,sixptcnctrdisc1z,sixptcnctrdisc1SPhi,sixptcnctrdisc1DPhi);
  sixptcnctrdisc1_log = new G4LogicalVolume(sixptcnctrdisc1,sixptcnctrdisc1Material,"sixptcnctrdisc1_log",0,0,0);
  //  sixptcnctrdisc1_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*sixptcnctrdisc1Pos),sixptcnctrdisc1_log,"sixptcnctrdisc1",expHall_log,false,0);
  LHTarget->AddPlacedVolume(sixptcnctrdisc1_log, *sixptcnctrdisc1Pos, &Cnctr1Rot);

  subcross2Rot = G4RotationMatrix::IDENTITY;
  subcross2Rot.rotateY(90.0*deg);
    
  sixptcnctr1 = new G4Tubs("sixptcnctr",sixptcnctrrmin,sixptcnctrrmax,sixptcnctrz,sixptcnctrSPhi,sixptcnctrDPhi);
  subcross1 = new G4Tubs("subcross1",subcross1rmin,subcross1rmax,subcross1z,subcross1SPhi,subcross1DPhi);
  subcross2 = new G4Tubs("subcross2",subcross2rmin,subcross2rmax,subcross2z,subcross2SPhi,subcross2DPhi);
  sixptcnctr2 = new G4SubtractionSolid("sixptcnctr2",sixptcnctr1,subcross1,G4Transform3D(NoRot,G4ThreeVector(0.,subcross1y,0.)));
  sixptcnctr3 = new G4SubtractionSolid("sixptcnctr3",sixptcnctr2,subcross2,G4Transform3D(subcross2Rot,G4ThreeVector(0.,subcross2y,0.)));
  sixptcnctr_log = new G4LogicalVolume(sixptcnctr3,sixptcnctrMaterial,"sixptcnctr_log",0,0,0);
  //  sixptcnctr_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*sixptcnctrPos),sixptcnctr_log,"sixptcnctr",expHall_log,false,0);
  LHTarget->AddPlacedVolume(sixptcnctr_log, *sixptcnctrPos, &Cnctr1Rot);

  cross1 = new G4Tubs("cross1",cross1rmin,cross1rmax,cross1z,cross1SPhi,cross1DPhi);
  cross1_ = new G4SubtractionSolid("cross1_",cross1,subcross2,G4Transform3D(subcross2Rot,G4ThreeVector(0.,subcross2y,0.)));
  cross1_log = new G4LogicalVolume(cross1_,cross1Material,"cross1_log",0,0,0);
  //  cross1_phys = new G4PVPlacement(G4Transform3D(NoRot,*cross1Pos),cross1_log,"cross1",expHall_log,false,0);
  LHTarget->AddPlacedVolume(cross1_log, *cross1Pos, &NoRot);

  cross2 = new G4Tubs("cross2",cross2rmin,cross2rmax,cross2z,cross2SPhi,cross2DPhi);
  cross2_log = new G4LogicalVolume(cross2,cross2Material,"cross2_log",0,0,0);
  //  cross2_phys = new G4PVPlacement(G4Transform3D(subcross2Rot,*cross2Pos),cross2_log,"cross2",expHall_log,false,0);
  LHTarget->AddPlacedVolume(cross2_log, *cross2Pos, &subcross2Rot);

  sixptcnctrdisc2 = new G4Tubs("sixptcnctrdisc2",sixptcnctrdisc2rmin,sixptcnctrdisc2rmax,sixptcnctrdisc2z,sixptcnctrdisc2SPhi,sixptcnctrdisc2DPhi);
  sixptcnctrdisc2_log = new G4LogicalVolume(sixptcnctrdisc2,sixptcnctrdisc2Material,"sixptcnctrdisc2_log",0,0,0);
  //  sixptcnctrdisc2_phys = new G4PVPlacement(G4Transform3D(Cnctr1Rot,*sixptcnctrdisc2Pos),sixptcnctrdisc2_log,"sixptcnctrdisc2",expHall_log,false,0);
  LHTarget->AddPlacedVolume(sixptcnctrdisc2_log, *sixptcnctrdisc2Pos, &Cnctr1Rot);
  
  LHTargetRot = G4RotationMatrix::IDENTITY;
  LHTargetRot.rotateZ(targetAngle);

  LHTarget->MakeImprint(expHall_log, *Pos, &LHTargetRot);

  G4Colour red (1.0,0.0,0.0); 
  G4VisAttributes* Vis_6 = new G4VisAttributes(red);
  Vis_6->SetVisibility(true);
  Vis_6->SetForceSolid(true);

  Target_log->SetVisAttributes(Vis_6);
 
  return Target_phys;
}

//-----------------------------------------------------------------------------
void Target::setX(G4double X)
{
  G4cout<<"----> Warning: User can't set Target X position of LH Target. Ignoring Target::setX("<< X << ")." <<G4endl;
}
//-----------------------------------------------------------------------------
void Target::setY(G4double Y)
{
   G4cout<<"----> Warning: User can't set Target Y position of LH Target. Ignoring Target::setY("<< Y << ")."<<G4endl;
}
//-----------------------------------------------------------------------------
void Target::setZ(G4double Z)
{
   G4cout<<"----> Warning: User can't set thickness of LH Target. Ignoring Target::setZ("<< Z << ")."<<G4endl;
}
//-----------------------------------------------------------------------------
void Target::setNStep(G4int n)
{
   NStep=n;
   target_limits->SetMaxAllowedStep(Target_thickness/NStep);
   Target_log->SetUserLimits(target_limits);
   G4cout<<"----> Number of simulation steps in the target is set to "<<NStep<<G4endl;
}

//-----------------------------------------------------------------------------
void Target::Report()
{
  G4cout<<"----> Target material set to  "<<Target_log->GetMaterial()->GetName()<< G4endl;   
  G4cout<<"----> Target density:         "<<Target_log->GetMaterial()->GetDensity()<< G4endl;   
  G4cout<<"----> Target radius set to "<<G4BestUnit(Target_rmax,"Length")<< G4endl;
  G4cout<<"----> Number of simulation steps in the target is set to "<<NStep<<G4endl;
}
//---------------------------------------------------------------------
void Target::setMaterial(G4String materialName)
{
   G4cout<<"----> Warning: User can't set the material of LH Target. Ignoring Target::setMaterial("<< materialName << ")."<<G4endl;
}
//-------------------------------------------------------------------
void Target::setTargetReactionDepth(G4double depth)
{
  //  G4cout<<"\n----> The depth is "<<G4BestUnit(depth,"Length")<< G4endl;;
  target_limits->SetUserMinRange(depth);
}
//-----------------------------------------------------------------------------
void Target::SetPositionZ(G4double d)
{
   G4cout<<"----> Warning: User can't set Z position of LH Target. Ignoring Target::setPositionZ("<< d << ")."<<G4endl;
}
//---------------------------------------------------------------------
void Target::ScaleDensity(G4double scale)
{
  // search the material by its name 
  G4String name=TargetMaterial->GetName();
  G4double Z=TargetMaterial->GetZ();
  G4double A=TargetMaterial->GetA();
  G4double density=TargetMaterial->GetDensity();
  density*=scale;
  TargetMaterial=new G4Material(name, Z,A,density);
  Target_log->SetMaterial(TargetMaterial);
  G4cout<<"----> Target material set to     "<<Target_log->GetMaterial()->GetName()<< G4endl;  
  G4cout<<"----> Target Z set to            "<<Target_log->GetMaterial()->GetZ()<< G4endl;  
  G4cout<<"----> Target mole mass set to       "<<Target_log->GetMaterial()->GetA()/g*mole<<" g/mole"<< G4endl;  
  G4cout<<"----> Target density set to         "<<Target_log->GetMaterial()->GetDensity()/g*cm3<<" g/cm3"<< G4endl;     
             
}
//---------------------------------------------------------------------
void Target::setSourceFrame(G4String sF)
{
  sourceFrame = sF;

  if(sourceFrame == "eu152_Z2707"){

    frameMaterial = materials->FindMaterial("Al");
    frameThickness = 2.9*mm;
    frameInnerRadius = 3.8*cm/2.0;
    frameOuterRadius = 5.4*cm/2.0;
    tapeMaterial = materials->FindMaterial("G4_POLYETHYLENE");
    tapeThickness = 0.012*cm;

    tape_r = frameInnerRadius - .2*cm; //in order to remove overlap

    euFrame = new G4Tubs("euFrame",frameInnerRadius,frameOuterRadius,frameThickness/2.,0.,360.*deg);
    euFrame_log = new G4LogicalVolume(euFrame,frameMaterial,"euFrame_log",0,0,0);
    euFrame_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),euFrame_log,"euFrame",expHall_log,false,0);

    euTape = new G4Tubs("euTape",0.,tape_r,tapeThickness/2.,0.,360.*deg);
    euTape_log = new G4LogicalVolume(euTape,tapeMaterial,"euTape_log",0,0,0);
    euTape_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),euTape_log,"euTape",expHall_log,false,0);

  } else if(sourceFrame == "cs137_E2879"){

    frameMaterial = materials->FindMaterial("Al");
    frameThickness = 0.7*mm;
    frameInnerRadius = 2.54*13./16.*cm/2.0;
    frameOuterRadius = 2.54*cm/2.0;
    tapeMaterial = materials->FindMaterial("G4_POLYETHYLENE");
    tapeThickness = 0.016*cm;

    csFrame = new G4Tubs("csFrame",frameInnerRadius,frameOuterRadius,frameThickness/2.,0.,360.*deg);
    csFrame_log = new G4LogicalVolume(csFrame,frameMaterial,"csFrame_log",0,0,0);
    csFrame_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),csFrame_log,"csFrame",expHall_log,false,0);

    G4ThreeVector ringPos(0.,0.,frameOuterRadius-frameInnerRadius+frameThickness/2.0);

    csRing = new G4Tubs("csRing",frameOuterRadius-frameThickness,frameOuterRadius,frameOuterRadius-frameInnerRadius,0.,360.*deg);
    csRing_log = new G4LogicalVolume(csRing,frameMaterial,"csRing_log",0,0,0);
    csRing_phys = new G4PVPlacement(G4Transform3D(NoRot,ringPos),csRing_log,"csRing",expHall_log,false,0);

    // Let's assume Mylar = Kapton = polyethylene is good enough.
    csTape = new G4Tubs("csTape",0.,frameInnerRadius,tapeThickness/2.,0.,360.*deg);
    csTape_log = new G4LogicalVolume(csTape,tapeMaterial,"csTape_log",0,0,0);
    csTape_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),csTape_log,"csTape",expHall_log,false,0);

  }

  G4cout<<"----> Source frame is set to "<<sourceFrame<< G4endl;                 
}
//-------------------------------------------------------------------
void Target::setSled()
{

  G4cout<<"----> Warning: target sled specified with LH target. Proceeding with no sled. "<< G4endl;                 

}
//-------------------------------------------------------------------
#endif
