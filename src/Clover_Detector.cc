#ifdef SCANNING
#include "Clover_Detector.hh"

Clover_Detector::Clover_Detector(G4LogicalVolume* experimentalHall_log,
				 Materials* mat, G4String orient)
{
  orientation = orient;

  materials=mat;
  expHall_log=experimentalHall_log;

  HpGe = materials->FindMaterial("HpGe");
  Al = materials->FindMaterial("Al");
  Cu = materials->FindMaterial("Cu");

  Length      = 8.0*cm;  // crystal length
  Radius      = 2.5*cm;  // crystal radius
  boxlength   = 4.5*cm;
  torusradius = 0.5*cm;  // "bevel" radius of front-face edges
  covergap    = 0.5*cm;  // front face of cover to front faces of crystals

  startAngle    = 0.*deg;
  spanningAngle = 360.*deg;

  Pos0.setX(0);
  Pos0.setY(0);
  Pos0.setZ(0);

  thetad = 90.*deg;
  phid = 90.*deg;

  Rot0=G4RotationMatrix::IDENTITY;
  Rot0.rotateX(180.*deg);
  Rot0.rotateY(90.*deg+thetad);

  LeafShift = 2.23*cm; // x and y offset relative to central axis

  // Lower left (facing clover)
  Leaf0Shift.setX(-LeafShift);
  Leaf0Shift.setY(-LeafShift);
  Leaf0Shift.setZ((Length + torusradius)/2. + covergap);
  Leaf0Pos = Pos0 + Leaf0Shift;

  // Lower right (facing clover)
  Leaf1Shift.setX(LeafShift);
  Leaf1Shift.setY(-LeafShift);
  Leaf1Shift.setZ((Length + torusradius)/2. + covergap);
  Leaf1Pos = Pos0 + Leaf1Shift;

  // Upper right (facing clover)
  Leaf2Shift.setX(LeafShift);
  Leaf2Shift.setY(LeafShift);
  Leaf2Shift.setZ((Length + torusradius)/2. + covergap);
  Leaf2Pos = Pos0 + Leaf2Shift;

  // Upper left (facing clover)
  Leaf3Shift.setX(-LeafShift);
  Leaf3Shift.setY(LeafShift);
  Leaf3Shift.setZ((Length + torusradius)/2. + covergap);
  Leaf3Pos = Pos0 + Leaf3Shift;

  CCoffset = 0.06*cm; // central contact x and y offset relative to the box
  CCradius = .6*cm; // central contact radius

  wallrot=G4RotationMatrix::IDENTITY;
  
  wallZoffset = 10.0*cm;

  coverlength    = 24.0*cm;
  coverwidth     = 10.1*cm;
  coverthickness = .2*cm;
  covershift.setZ(coverlength/2.);
  coverpos = Pos0 + covershift;

  cornerRadius = 1.55*cm;
  corneroffset = coverwidth/2. - cornerRadius;

  Cuboxlength = 4.00*cm;

  cornershift.setX(-corneroffset);
  cornershift.setY(corneroffset);
  cornershift.setZ(wallZoffset);
  cornerpos = Pos0 + cornershift;

  corner1shift.setX(-corneroffset);
  corner1shift.setY(-corneroffset);
  corner1shift.setZ(wallZoffset);
  corner1pos = Pos0 + corner1shift;

  corner2shift.setX(corneroffset);
  corner2shift.setY(-corneroffset);
  corner2shift.setZ(wallZoffset);
  corner2pos = Pos0 + corner2shift;

  corner3shift.setX(corneroffset);
  corner3shift.setY(corneroffset);
  corner3shift.setZ(wallZoffset);
  corner3pos = Pos0 + corner3shift;

  Cuboxshift.setZ(Length + covergap + Cuboxlength/2.);
  Cuboxpos = Pos0 + Cuboxshift;

  // Final Clover placement
  if(orientation == "right")
    DetPos.setX(  -72.06*mm );
  else if(orientation == "left")
    DetPos.setX(   72.06*mm );
  DetPos.setY(  343.88*mm ); // Shield y position
  DetPos.setZ( -155.44*mm );
  
  DetRot=G4RotationMatrix::IDENTITY;
  if(orientation == "right")
    DetRot.rotateY(200.*deg);
  else if(orientation == "left")
    DetRot.rotateY(160.*deg);
}

Clover_Detector::~Clover_Detector()
{
}

G4VPhysicalVolume* Clover_Detector::Construct()
{

  // Leaf

  detector = new G4Tubs("detector", 0, Radius, (Length-torusradius)/2., 
			0., 360.*deg);

  CCsub = new G4Tubs("CCsub", 0, CCradius, Length/2., 0., 360.*deg);

  box = new G4Box("box", boxlength/2., boxlength/2., 9.0*cm);

  // creating curved bevel on the front edge

  torus = new G4Torus("torus", 0, torusradius, Radius - torusradius,
		      0., 360.*deg);

  torusbox = new G4Box("torusbox", Radius + torusradius*3., 
		       Radius+torusradius*3., 
		       torusradius);
 
  torustube = new G4Tubs("torustube", 0, Radius-torusradius, torusradius, 
			 0., 360.*deg);

  torus1 = new G4UnionSolid("torus2", torustube, torus,
			    G4Transform3D(G4RotationMatrix(),G4ThreeVector()));

  torus2 = 
    new G4SubtractionSolid("torus1", torus1, torusbox, 
			   G4Transform3D(G4RotationMatrix(),
					 G4ThreeVector(0.,0.,torusradius)));

  //Actually building detector

  bevel = 
    new G4UnionSolid("bevel", detector, torus2,
		     G4Transform3D(G4RotationMatrix(),
				   G4ThreeVector(0,0,
						 -(Length-torusradius)/2.)));

  subtract = 
    new G4SubtractionSolid("subtraction", bevel, CCsub,
			   G4Transform3D(G4RotationMatrix(),
					 G4ThreeVector(0.*cm,0.*cm,1.*cm)));

  intersect = new G4IntersectionSolid ("intersect", subtract, box, 
				       G4Transform3D(G4RotationMatrix(),
						     G4ThreeVector(CCoffset, 
								   CCoffset, 
								   0.)));

  detector_log = new G4LogicalVolume(intersect, HpGe, "Leaf_log", 0, 0, 0);

  // Cryosatat

  boxout = new G4Box("box", coverwidth/2., coverwidth/2., coverlength/2.);

  boxin = new G4Box("boxin", coverwidth/2. - coverthickness, 
		    coverwidth/2. - coverthickness, coverlength/2.);

  cover = 
    new G4SubtractionSolid("cover", boxout, boxin, 
			   G4Transform3D(G4RotationMatrix(),
					 G4ThreeVector(0,0,coverthickness)));

  // Make space for the curved corners.

  cornerCut = new G4Tubs("cornerCut", cornerRadius - coverthickness,
			 cornerRadius*4., coverlength, 0*deg, 90.*deg);

  coversub = new G4SubtractionSolid("cover1", cover, cornerCut,
				  G4Transform3D(G4RotationMatrix(),
						G4ThreeVector(corneroffset, 
							      corneroffset, 
							      0.*cm)));

  coversub = 
    new G4SubtractionSolid("cover2", coversub, cornerCut,
			   G4Transform3D(G4RotationMatrix(0., 0., 90.*deg),
					 G4ThreeVector(corneroffset, 
						       -corneroffset,
						       0.)));

  coversub = 
    new G4SubtractionSolid("cover3", coversub, cornerCut,
			   G4Transform3D(G4RotationMatrix(0., 0., 180.*deg),
					 G4ThreeVector(-corneroffset,
						       -corneroffset,
						       0.)));

  coversub = 
    new G4SubtractionSolid("cover4", coversub, cornerCut,
			   G4Transform3D(G4RotationMatrix(0., 0., 270.*deg),
					 G4ThreeVector(-corneroffset,
						       corneroffset,
						       0.)));

  // Add the corners.

  corner = new G4Tubs("corner", cornerRadius - coverthickness,
		      cornerRadius, coverlength/2., 0., 90.*deg);

  coveru = new G4UnionSolid("cover5", coversub, corner, 
			    G4Transform3D(G4RotationMatrix(),
					  G4ThreeVector(corneroffset,
							corneroffset,
							0.)));

  coveru = new G4UnionSolid("cover6", coveru, corner, 
			    G4Transform3D(G4RotationMatrix(0., 0., 90.*deg),
					  G4ThreeVector(corneroffset,
							-corneroffset,
							0.)));

  coveru = new G4UnionSolid("cover7", coveru, corner,
			    G4Transform3D(G4RotationMatrix(0., 0., 180.*deg),
					  G4ThreeVector(-corneroffset,
							-corneroffset,
							0.)));
  
  coveru = new G4UnionSolid("cover8", coveru, corner,
			    G4Transform3D(G4RotationMatrix(0., 0., 270.*deg),
					  G4ThreeVector(-corneroffset,
							corneroffset,
							0.)));

  cover_log = new G4LogicalVolume(coveru, Al, "cover_log", 0, 0, 0);

  // Copper Backing

  Cubox = new G4Box("copperbox", boxlength + .03*cm, boxlength + .03*cm, 
		    Cuboxlength/2.);

  CuboxCut = new G4Tubs("CuboxCut", Radius, Radius*1.5, Cuboxlength + 1.*cm, 
			0., 90.*deg); 

  CuboxCut1 = 
    new G4SubtractionSolid("CuboxCut1", Cubox, CuboxCut,
			   G4Transform3D(G4RotationMatrix(0., 0., 0.),
					 G4ThreeVector(LeafShift,
						       LeafShift,
						       0.)));
  
  CuboxCut2 = 
    new G4SubtractionSolid("CuboxCut2", CuboxCut1, CuboxCut,
			   G4Transform3D(G4RotationMatrix(0., 0., 90.*deg),
					 G4ThreeVector(LeafShift,
						       -LeafShift,
						       0.)));
 
  CuboxCut3 = 
    new G4SubtractionSolid("CuboxCut3", CuboxCut2, CuboxCut,
			   G4Transform3D(G4RotationMatrix(0., 0., 180.*deg),
					 G4ThreeVector(-LeafShift,
						       -LeafShift,
						       0.)));

  CuboxCut4 = 
    new G4SubtractionSolid("CuboxCut4", CuboxCut3, CuboxCut, 
			   G4Transform3D(G4RotationMatrix(0., 0., 270.*deg),
					 G4ThreeVector(-LeafShift,
						       LeafShift,
						       0.)));

  Cubox_log = new G4LogicalVolume(CuboxCut4,Cu,"Cubox_log",0,0,0);

  assemblyclover = new G4AssemblyVolume();

  assemblyclover->AddPlacedVolume(detector_log,Leaf0Pos,&Rot0);
 
  Rot0.rotateZ(90.*deg);

  assemblyclover->AddPlacedVolume(detector_log,Leaf1Pos,&Rot0);

  Rot0.rotateZ(90.*deg);
  
  assemblyclover->AddPlacedVolume(detector_log,Leaf2Pos,&Rot0);

  Rot0.rotateZ(90.*deg);

  assemblyclover->AddPlacedVolume(detector_log,Leaf3Pos,&Rot0);

  assemblyclover->AddPlacedVolume(cover_log,coverpos,&wallrot);

  assemblyclover->AddPlacedVolume(Cubox_log,Cuboxpos,&wallrot);

  assemblyclover->MakeImprint(expHall_log, DetPos, &DetRot);

  //Visualization Attributes

  //   Clover Crystal
  G4Colour dgreen (0.0,0.75, 0.0, 1.0); 
  G4VisAttributes* Vis_1 = new G4VisAttributes(dgreen);
  Vis_1->SetVisibility(true);
  Vis_1->SetForceSolid(true);

  //   Can
  G4Colour transGrey (0.8, 0.8, 0.8, 0.2);
  G4VisAttributes* Vis_2 = new G4VisAttributes(transGrey);
  Vis_2->SetVisibility(true);
  Vis_2->SetForceSolid(false);
  Vis_2->SetForceWireframe(true);

  detector_log->SetVisAttributes(Vis_1);
  cover_log->SetVisAttributes(Vis_2);
  Cubox_log->SetVisAttributes(Vis_2);
 
  return detector_phys;
}
//---------------------------------------------------------------------
void Clover_Detector::MakeSensitive(TrackerGammaSD* TrackerGamma)
{
  detector_log->SetSensitiveDetector(TrackerGamma);
}
//---------------------------------------------------------------------
void Clover_Detector::setX(G4double x)
{
  DetPos.setX(x);
  //  detector_phys->SetTranslation(DetPos);
}
//---------------------------------------------------------------------
void Clover_Detector::setY(G4double y)
{
  DetPos.setY(y);
  //  detector_phys->SetTranslation(DetPos);
}
//---------------------------------------------------------------------
void Clover_Detector::setZ(G4double z)
{
  DetPos.setZ(z);
  //  detector_phys->SetTranslation(DetPos);
}
//---------------------------------------------------------------------
#endif
