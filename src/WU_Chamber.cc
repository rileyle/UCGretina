#ifndef LHTARGET
#include "WU_Chamber.hh"

WU_Chamber::WU_Chamber(Materials* mat)
{
  materials=mat;
  BTrmin = 3.0*2.54*cm - 0.058*2.54*cm;
  BTrmax = 3.0*2.54*cm;
  BTDz=44.*cm; //LR (approx target to gate valve flange)
  BTSPhi=0.*deg;
  BTDPhi=360*deg; 
  pRmin=171*mm;
  pRmax=173*mm;
  pSPhi=-90*deg;
  pDPhi=180*deg;
  pSTheta=0.*deg;
  pDTheta=180*deg;
  exBox=10.*mm;
  smallCap=52.925*mm;
  cRmin = 156*mm;
  cRmax = 158*mm;
  boxSide = 31.75*mm;
  boxSideLarge = 34.75*mm;
  sine1=141.3133925*mm;
  cosine1=72.67195405*mm;
  apertureA=25.4*mm;
  apertureB=25.4*mm;
  apertureC=25.4*mm;
  RingDistance = cRmax + 106.92*mm;
  LTRadius1=38.2/2*mm;
  LTRadius2=57.15/2*mm;
  LTRadius3=53.15/2*mm;
  LTLength1=28.5*mm;
  LTLength2=120*mm + 84.33*mm;
  LTDistance= RingDistance;
  Pos0 = new G4ThreeVector(0.,0.,0.);
  BeamTubeMaterial = materials->FindMaterial("Al");
  Brass = materials->FindMaterial("Cu");
  Glass = materials->FindMaterial("Si");
}

WU_Chamber::~WU_Chamber()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* WU_Chamber::Construct(G4LogicalVolume* experimentalHall_log)
{
  expHall_log=experimentalHall_log;
  //Building the hemispheres
  Hemisphere = new G4Sphere("Hemisphere",pRmin,pRmax,pSPhi,pDPhi,pSTheta,pDTheta);

  G4Tubs* Cylinder = new G4Tubs( "Cylinder", 0., 1.1*pRmax, smallCap, 0.*deg, 360.*deg);

  G4RotationMatrix Rot = G4RotationMatrix::IDENTITY;
  Rot.rotateY(90.*deg);

  G4SubtractionSolid* Cap = new G4SubtractionSolid("Cap", Hemisphere, Cylinder, G4Transform3D(Rot,G4ThreeVector(0.0,0.0,0.0)));

  Cap_log = new G4LogicalVolume(Cap,BeamTubeMaterial,"Cap_log",0,0,0);

  Cap_phys1 = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
			   Cap_log,"ChamberCap1",expHall_log,false,0);
  
  G4RotationMatrix myNewRot = G4RotationMatrix::IDENTITY;
  myNewRot.rotateZ(180.*deg);

  Cap_phys2 = new G4PVPlacement(G4Transform3D(myNewRot,*Pos0),
			    Cap_log,"ChamberCap2",expHall_log,false,0);

////////////////////////////////////////////////////////////////////////////////

  CylinderFill = new G4Tubs( "CylinderFill", 0.*mm, cRmax, smallCap, 0.*deg, 360.*deg); //building cylinder in between hemispheres

  G4Box* Box = new G4Box( "Box", boxSide, boxSide, exBox); //builds raised boxes on outer surface of cylinder (same code for all 6 boxes)

  G4RotationMatrix BoxRot = G4RotationMatrix::IDENTITY;
  BoxRot.rotateY(90.*deg);

  G4UnionSolid* CylBox1 = new G4UnionSolid("CylBox1", CylinderFill, Box, G4Transform3D(BoxRot,G4ThreeVector(cRmax,0.0,0.0)));

  CylBox1_log = new G4LogicalVolume(CylBox1,BeamTubeMaterial,"CylBox1_log",0,0,0);

  G4Box* Box2 = new G4Box( "Box2", boxSide, boxSide, exBox);  

  G4RotationMatrix BoxRot2 = G4RotationMatrix::IDENTITY;
  BoxRot2.rotateX(90.*deg);
  BoxRot2.rotateZ(153.43*deg);

  G4UnionSolid* CylBox2 = new G4UnionSolid("CylBox2", CylBox1, Box2, G4Transform3D(BoxRot2,G4ThreeVector(cosine1,sine1,0.0)));

  CylBox2_log = new G4LogicalVolume(CylBox2,BeamTubeMaterial,"CylBox2_log",0,0,0);

  G4Box* Box3 = new G4Box( "Box3", boxSide, boxSide, exBox);  

  G4RotationMatrix BoxRot3 = G4RotationMatrix::IDENTITY;
  BoxRot3.rotateX(90.*deg);

  G4UnionSolid* CylBox3 = new G4UnionSolid("CylBox3", CylBox2, Box3, G4Transform3D(BoxRot3,G4ThreeVector(0.0,cRmax,0.0)));

  CylBox3_log = new G4LogicalVolume(CylBox3,BeamTubeMaterial,"CylBox3_log",0,0,0);

  G4Box* Box4 = new G4Box( "Box4", boxSideLarge, boxSideLarge, exBox);  

  G4RotationMatrix BoxRot4 = G4RotationMatrix::IDENTITY;
  BoxRot4.rotateX(90.*deg);
  BoxRot4.rotateZ(63.43*deg);

  G4UnionSolid* CylBox4 = new G4UnionSolid("CylBox4", CylBox3, Box4, G4Transform3D(BoxRot4,G4ThreeVector(-sine1,cosine1,0.0)));

  CylBox4_log = new G4LogicalVolume(CylBox4,BeamTubeMaterial,"CylBox4_log",0,0,0);

  G4Box* Box5 = new G4Box( "Box5", boxSide, boxSide, exBox);

  G4UnionSolid* CylBox5 = new G4UnionSolid("CylBox5", CylBox4, Box5, G4Transform3D(BoxRot,G4ThreeVector(-cRmax,0.0,0.0)));

  CylBox5_log = new G4LogicalVolume(CylBox5,BeamTubeMaterial,"CylBox5_log",0,0,0);

  G4Box* Box6 = new G4Box( "Box6", boxSide, boxSide, exBox);

  G4UnionSolid* CylBox6 = new G4UnionSolid("CylBox6", CylBox5, Box6, G4Transform3D(BoxRot2,G4ThreeVector(-cosine1,-sine1,0.0)));

  CylBox6_log = new G4LogicalVolume(CylBox6,BeamTubeMaterial,"CylBox6_log",0,0,0);

////////////////////////////////////////////////////////////////////////////////

  G4Tubs* Tube1 = new G4Tubs( "Tube1", 0., apertureA, 1.1*pRmax, 0.*deg, 360.*deg); //Subtraction solid through the first and fifth raised boxes to accomodate the incoming beam tube

  G4RotationMatrix TubeRot1 = G4RotationMatrix::IDENTITY;
  TubeRot1.rotateY(90.*deg);

  G4SubtractionSolid* BeamHole1 = new G4SubtractionSolid("BeamHole1", CylBox6, Tube1, G4Transform3D(TubeRot1,G4ThreeVector(0.0,0.0,0.0)));

  BeamHole1_log = new G4LogicalVolume(BeamHole1,BeamTubeMaterial,"BeamHole1_log",0,0,0);

  G4Tubs* Tube2 = new G4Tubs( "Tube2", 0., apertureA, 1.1*pRmax, 0.*deg, 360.*deg); //Subtraction solid through second and sixth raised boxes for target cable or ladder port 

  G4RotationMatrix TubeRot2 = G4RotationMatrix::IDENTITY;
  TubeRot2.rotateX(90.*deg);
  TubeRot2.rotateZ(153.43*deg);

  G4SubtractionSolid* BeamHole2 = new G4SubtractionSolid("BeamHole2", BeamHole1, Tube2, G4Transform3D(TubeRot2,G4ThreeVector(0.0,0.0,0.0)));

  BeamHole2_log = new G4LogicalVolume(BeamHole2,BeamTubeMaterial,"BeamHole2_log",0,0,0);

  G4Tubs* Tube3 = new G4Tubs( "Tube3", 0., apertureB, 0.5*pRmax, 0.*deg, 360.*deg); //Subtraction solid through third raised box to accomodate target port

  G4RotationMatrix TubeRot3 = G4RotationMatrix::IDENTITY;
  TubeRot3.rotateX(90.*deg);

  G4SubtractionSolid* BeamHole3 = new G4SubtractionSolid("BeamHole3", BeamHole2, Tube3, G4Transform3D(TubeRot3,G4ThreeVector(0.0,cRmax,0.0)));

  BeamHole3_log = new G4LogicalVolume(BeamHole3,BeamTubeMaterial,"BeamHole3_log",0,0,0);

  G4Tubs* Tube4 = new G4Tubs( "Tube4", 0., apertureC, 0.5*pRmax, 0.*deg, 360.*deg); //Subtraction solid through fourth raised box for target viewing port

  G4RotationMatrix TubeRot4 = G4RotationMatrix::IDENTITY;
  TubeRot4.rotateX(90.*deg);
  TubeRot4.rotateZ(63.43*deg);

  G4SubtractionSolid* BeamHole4 = new G4SubtractionSolid("BeamHole4", BeamHole3, Tube4, G4Transform3D(TubeRot4,G4ThreeVector(-sine1,cosine1,0.0)));

  BeamHole4_log = new G4LogicalVolume(BeamHole4,BeamTubeMaterial,"BeamHole4_log",0,0,0);

  G4Tubs* InnerCylinder = new G4Tubs( "InnerCylinder", 0.*mm, cRmin, 1.1*smallCap, 0.*deg, 360.*deg); //Subtraction solid to hollow out the cylinder in order to get it as a shell

  G4SubtractionSolid* InCyl = new G4SubtractionSolid("InCyl", BeamHole4, InnerCylinder, G4Transform3D(NoRot,G4ThreeVector(0.0,0.0,0.0)));

  InCyl_log = new G4LogicalVolume(InCyl,BeamTubeMaterial,"InCyl_log",0,0,0);

  CylBox_phys = new G4PVPlacement(G4Transform3D(Rot,*Pos0),
			    InCyl_log,"ChamberInCyl",expHall_log,false,0);

  
////////////////////////////////////////////////////////////////////////////////
//Connection Tube

  const G4double CzPlane[21] = {0., 8.89*mm, 8.89*mm, 12.55*mm, 12.55*mm, 13.29*mm, 13.29*mm, 14.14*mm, 14.14*mm, 14.48*mm, 84.33*mm, 84.33*mm, 86.36*mm, 86.36*mm, 104.18*mm, 104.18*mm, 113.68*mm, 113.68*mm, 118.45*mm, 118.45*mm, 123.2*mm};
  const G4double CrInner[21] = {44.75/2*mm, 44.75/2*mm, 44.75/2*mm, 35.475/2*mm, 35.475/2*mm, 33.60/2*mm, 33.60/2*mm, 33.60/2*mm, 33.60/2*mm, 33.60/2*mm, 33.60/2*mm, 33.60/2*mm, 33.60/2*mm, 33.60/2*mm, 44.45/2*mm, 44.45/2*mm, 44.45/2*mm, 44.45/2*mm, 41.4/2*mm, 41.4/2*mm, 46.94/2*mm};
  const G4double CrOuter[21] = {50.8/2*mm, 50.8/2*mm, 50.8/2*mm, 63.5/2*mm, 63.5/2*mm, 52.036/2*mm, 52.036/2*mm, 36.60/2*mm, 36.60/2*mm, 36.60/2*mm, 36.60/2*mm, 36.60/2*mm, 47.63/2*mm, 47.63/2*mm, 47.63/2*mm, 47.63/2*mm, 50.8/2*mm, 50.8/2*mm, 50.8/2*mm, 50.8/2*mm, 50.8/2*mm};

  Connector = new G4Polycone( "Connector", 0., 360.*deg, 
				21, CzPlane, CrInner,  CrOuter);  //Constructs the tube that will be connected to the first and fifth beam 									    ports

  G4RotationMatrix ConRot = G4RotationMatrix::IDENTITY;
  ConRot.rotateY(180.*deg);

  Connector_log = new G4LogicalVolume(Connector,BeamTubeMaterial,"Connnector_log",0,0,0);

  Connector_phys 
    = new G4PVPlacement(G4Transform3D(ConRot,
				      G4ThreeVector(0.0,0.0,-cRmax)),
			Connector_log,"ChamberConnector",expHall_log,
			false,0);

  Connector2_phys 
    = new G4PVPlacement(G4Transform3D(NoRot,
				      G4ThreeVector(0.0,0.0,cRmax)),
			Connector_log,"ChamberConnector2",expHall_log,
			false,0);


  const G4double RzPlane[8] = {0., 2.375*mm,  2.375*mm, 4.75*mm, 4.75*mm, 12.67*mm, 12.67*mm, 19.03*mm};

  const G4double RrInner[8] = {50.8/2*mm, 50.8/2*mm, 50.8/2*mm, 50.8/2*mm, 50.8/2*mm, 41.4/2*mm, 41.4/2*mm, 41.4/2*mm};

  const G4double RrOuter[8] = {50.8/2*mm, 53.975/2*mm, 53.975/2*mm, 57.15/2*mm, 57.15/2*mm, 57.15/2*mm, 57.15/2*mm, 54.61/2*mm};

  Ring = new G4Polycone( "Ring", 0., 360.*deg, 8,  RzPlane, RrInner,  RrOuter);  // Constructs the brass ring that will connect the tube we 											    just constructed to the larger tube carrying the beam

  Ring_log = new G4LogicalVolume(Ring, Brass,"Ring_log",0,0,0);

  Ring_phys 
    = new G4PVPlacement(G4Transform3D(ConRot,
				      G4ThreeVector(0.0,0.0,-RingDistance)),
			Ring_log,"ChamberRing",expHall_log,
			false,0);

  Ring2_phys 
    = new G4PVPlacement(G4Transform3D(NoRot,
				      G4ThreeVector(0.0,0.0,RingDistance)),
			Ring_log,"ChamberRing2",expHall_log,
			false,0);

  const G4double LzPlane[6] = {0., LTLength1, LTLength1, LTLength1, LTLength1, LTLength2};

  const G4double LrInner[6] = {LTRadius1, LTRadius1, LTRadius1, LTRadius1, LTRadius1, LTRadius3};

  const G4double LrOuter[6] = {41.4/2*mm, LTRadius2, LTRadius2, LTRadius3, LTRadius3, LTRadius2};

  LongTube = new G4Polycone( "LongTube", 0., 360.*deg, 6,  LzPlane, LrInner,  LrOuter); //Constructs long tube that will be carrying the beam 												  to the cylinder

  LongTube_log = new G4LogicalVolume(LongTube, BeamTubeMaterial,"LongTube_log",0,0,0);

  LongTube_phys 
    = new G4PVPlacement(G4Transform3D(ConRot,
				      G4ThreeVector(0.0,0.0,-LTDistance)),
			LongTube_log,"ChamberLongTube",expHall_log,
			false,0);

  LongTube2_phys 
    = new G4PVPlacement(G4Transform3D(NoRot,
				      G4ThreeVector(0.0,0.0,LTDistance)),
			LongTube_log,"ChamberLongTube2",expHall_log,
			false,0);

  Window = new G4Tubs ("Window", 0., apertureC, exBox, 0.*deg, 360.*deg); //Constructs the small glass window that will fit into the target 										    viewing port

  G4RotationMatrix WindowRot = G4RotationMatrix::IDENTITY;
  WindowRot.rotateX(153.43*deg);

  Window_log = new G4LogicalVolume(Window,Glass,"Window_log",0,0,0);

  Window_phys 
    = new G4PVPlacement(G4Transform3D(WindowRot,
				      G4ThreeVector(0.0,cosine1,sine1)),
			Window_log,"ChamberWindow",expHall_log,
			false,0);

  const G4double PPzPlane[10] = {0., 8.89*mm, 8.89*mm, 8.89*mm, 8.89*mm, 12.55*mm, 12.55*mm, 13.55*mm, 13.55*mm, 14.48*mm};

  const G4double PPrInner[10] = {33.75/2*mm, 33.75/2*mm, 33.75/2*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm, 0.*mm};

  const G4double PPrOuter[10] = {39.8/2*mm, 39.8/2*mm, 39.8/2*mm, 63.5/2*mm, 63.5/2*mm, 63.5/2*mm, 63.5/2*mm, 55.603/2*mm, 55.603/2*mm, 36.6/2*mm};

  PortPlug = new G4Polycone( "PortPlug",  0., 360.*deg, 10, PPzPlane, PPrInner, PPrOuter); //Constructs the port plug which will go into the  												     cable or target ladder ports (second and sixth 												     raised boxes) 

  G4RotationMatrix PPRot = G4RotationMatrix::IDENTITY;
  PPRot.rotateX(243.43*deg);

  G4Tubs* SidePort = new G4Tubs( "SidePort", 0., 5.08/2*mm, 2.79*mm, 0.*deg, 360.*deg); //Makes the small cylindrical divets on the outside of 												  the port plug as a subtraction solid

  G4SubtractionSolid* Port1 = new G4SubtractionSolid("Port1", PortPlug, SidePort, 
                                                     G4Transform3D(NoRot,G4ThreeVector(13.97*mm,0.0,14.48*mm)));

  G4SubtractionSolid* Port = new G4SubtractionSolid("Port", Port1, SidePort, G4Transform3D(NoRot,G4ThreeVector(-13.97*mm,0.0,14.48*mm)));

  Port_log= new G4LogicalVolume(Port,BeamTubeMaterial,"Port_log",0,0,0);

  PortPlug_phys 
    = new G4PVPlacement(G4Transform3D(PPRot,
				      G4ThreeVector(0.0,sine1,-cosine1)),
			Port_log,"ChamberPort",expHall_log,
			false,0);

  G4RotationMatrix PPRot2 = G4RotationMatrix::IDENTITY;
  PPRot2.rotateX(63.43*deg);

  PortPlug2_phys 
    = new G4PVPlacement(G4Transform3D(PPRot2,
				      G4ThreeVector(0.0,-sine1,cosine1)),
			Port_log,"ChamberPort2",expHall_log,
			false,0);

  const G4double TPzPlane[16] = {0., 7.1214*mm, 7.1214*mm, 7.1214*mm, 7.1214*mm, 8.901*mm, 8.901*mm, 32.381*mm, 32.381*mm, 33.865*mm, 33.865*mm, 43.479*mm, 43.479*mm, 43.479*mm, 43.479*mm, 47.039*mm};

  const G4double TPrInner[16] = {6.35/2*mm, 6.35/2*mm, 6.35/2*mm, 6.35/2*mm, 6.35/2*mm, 6.35/2*mm, 6.35/2*mm, 44.2126/2*mm, 44.2126/2*mm, 44.51/2*mm, 44.51/2*mm, 44.51/2*mm, 44.51/2*mm, 44.51/2*mm, 44.51/2*mm, 44.51/2*mm};

  const G4double TPrOuter[16] = {11.869/2*mm, 11.869/2*mm, 11.869/2*mm, 22.5514/2*mm, 22.5514/2*mm, 24.33177/2*mm, 24.33177/2*mm, 48.66355/2*mm, 48.66355/2*mm, 50.8/2*mm, 50.8/2*mm, 50.8/2*mm, 50.8/2*mm, 63.51/2*mm, 63.51/2*mm, 63.51/2*mm};

  TargetPort = new G4Polycone( "TargetPort", 0., 360.*deg, 16,  TPzPlane, TPrInner,  TPrOuter); //Constructs the target port that fits into 													  the topmost (third) raised box

  TargetPort_log = new G4LogicalVolume(TargetPort, BeamTubeMaterial,"TargetPort_log",0,0,0);

  G4RotationMatrix TPRot = G4RotationMatrix::IDENTITY;
  TPRot.rotateX(270*deg);

  TargetPort_phys 
    = new G4PVPlacement(G4Transform3D(TPRot,
				      G4ThreeVector(0.0,cRmax-32.381*mm,0.0)),
			TargetPort_log,"ChamberTargetPort",expHall_log,
			false,0);

  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0,0.3);
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(false);

  Cap_log->SetVisAttributes(Vis);

  return CylBox_phys;
}
//-----------------------------------------------------------------------------
#endif
