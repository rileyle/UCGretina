#include "Clover_Detector.hh"

Clover_Detector::Clover_Detector(G4LogicalVolume* experimentalHall_log,
				 Materials* mat, G4String orient)
{
  orientation = orient;

  materials=mat;
  expHall_log=experimentalHall_log;

  HpGe = materials->FindMaterial("Germanium");
  Al = materials->FindMaterial("Al");
  Cu = materials->FindMaterial("Cu");

  // Eurysys CLOVER 4X50X80 SEG2 manual p. 19
  Length      = 80*mm;   // crystal length
  Radius      = 25*mm;   // crystal radius
  boxlength   = 45*mm;   // box cutting planar boundaries on crystals
  torusradius =  5*mm;   // "bevel" radius of front-face edges
  covergap    =  5*mm;   // front face of cover to front faces of crystals
  LeafShift   = 22.3*mm; // x and y crystal offset relative to central axis
  CCoffset    =  0.5*mm; // central contact x and y offset relative to the box

  CCradius    = 5*mm;    // central contact radius (based on GRETINA crystals)
  CCdepth     = 50*mm;   // David Radford/Mitch Allmond best guess

  coaxialDLThickness = 0.5*mm;
  outerDLThickness = 0.5*mm;
  
  startAngle    = 0.*deg;
  spanningAngle = 360.*deg;

  Pos0.setX(0);
  Pos0.setY(0);
  Pos0.setZ(0);

  thetad = 90.*deg;
  phid = 90.*deg;

  DetTheta = 0.; // For non scanning-table placement
  DetPhi = 0.;   // For non scanning-table placement

  DetCode = 0;   // For non scanning-table placement
  
  Rot0=G4RotationMatrix::IDENTITY;
  Rot0.rotateX(180.*deg);
  Rot0.rotateY(90.*deg+thetad);

  // Upper left (facing clover)
  Leaf0Shift.setX(LeafShift);
  Leaf0Shift.setY(LeafShift);
  Leaf0Shift.setZ((Length + torusradius)/2. + covergap);
  Leaf0Pos = Pos0 + Leaf0Shift;

  // Upper right (facing clover)
  Leaf1Shift.setX(-LeafShift);
  Leaf1Shift.setY(LeafShift);
  Leaf1Shift.setZ((Length + torusradius)/2. + covergap);
  Leaf1Pos = Pos0 + Leaf1Shift;

  // Lower right (facing clover)
  Leaf2Shift.setX(-LeafShift);
  Leaf2Shift.setY(-LeafShift);
  Leaf2Shift.setZ((Length + torusradius)/2. + covergap);
  Leaf2Pos = Pos0 + Leaf2Shift;

  // Lower left (facing clover)
  Leaf3Shift.setX( LeafShift);
  Leaf3Shift.setY(-LeafShift);
  Leaf3Shift.setZ((Length + torusradius)/2. + covergap);
  Leaf3Pos = Pos0 + Leaf3Shift;

  wallrot=G4RotationMatrix::IDENTITY;
  
  wallZoffset = 100*mm;

  coverlength    = 240*mm;
  coverwidth     = 101*mm;
  coverthickness = 1.5*mm;
  covershift.setZ(coverlength/2.);
  coverpos = Pos0 + covershift;

  cornerRadius = 15.5*mm;
  corneroffset = coverwidth/2. - cornerRadius;

  Cuboxlength = 40.0*mm;

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

  // Final Clover placement (scanning table)
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

void Clover_Detector::Construct()
{

  if(orientation = "FDS"){
    DetPos.rotateZ(DetPsi);
    DetPos.rotateY(DetTheta);
    DetPos.rotateZ(DetPhi);
    DetRot.rotateY(DetPos.getTheta());
    DetRot.rotateZ(DetPos.getPhi());
  }

  //Visualization Attributes ===================================================

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
  
  // Single clover leaf ========================================================

  G4IntersectionSolid* intersect = Bulletized_Boxed_Cylinder(0);
  
  // ... central contact
  CCsub = new G4Tubs("CCsub", 0, CCradius, Length/2., 0., 360.*deg);

  G4SubtractionSolid* leaf =
    new G4SubtractionSolid("subtraction", intersect, CCsub,
			   G4Transform3D(G4RotationMatrix(),
					 G4ThreeVector(0., 0,
						       Length-CCdepth)));

  leaf_log = new G4LogicalVolume(leaf, HpGe, "Leaf_log", 0, 0, 0);

  leaf_log->SetVisAttributes(Vis_1);
  
  // G4cout << "\n  Total HPGe volume = "
  // 	 <<std::fixed<<std::setprecision(3)<<std::setw(7)<<std::right
  // 	 << leaf->GetCubicVolume()/cm3*4.
  // 	 << " cm3" << G4endl;
  
  // Outer dead layer (including back DL for now)

  if(outerDLThickness > 0){
    G4IntersectionSolid* outerDLcutout =
      Bulletized_Boxed_Cylinder(outerDLThickness);

    G4SubtractionSolid* outerDL =
      new G4SubtractionSolid("subtraction", intersect, outerDLcutout,
			     G4Transform3D(G4RotationMatrix(),
					   G4ThreeVector(0., 0, 0.)));
    outerDL =
      new G4SubtractionSolid("outerDL", outerDL, CCsub,
			     G4Transform3D(G4RotationMatrix(),
					   G4ThreeVector(0., 0,
							 Length-CCdepth)));

    G4LogicalVolume* outerDL_log =
      new G4LogicalVolume( outerDL, HpGe, "outerDL_log", 0, 0, 0 );

    // Place outer dead layer as a daugter of the leaf logical volume
    new G4PVPlacement( 0, G4ThreeVector(), outerDL_log, "outerDeadLayer",
		       leaf_log, false, 0);

    outerDL_log->SetVisAttributes(Vis_2);
  }
  
  // Coaxial dead layer
  if(coaxialDLThickness > 0){
    G4double zSlice[4] = {Length/2.-CCdepth-coaxialDLThickness,
			  Length/2.-CCdepth,
			  Length/2.-CCdepth,
			  Length/2.-outerDLThickness};
    G4double rIn[4]    = {0., 0., CCradius, CCradius};
    G4double rOut[4]   = {CCradius+coaxialDLThickness,
			  CCradius+coaxialDLThickness,
			  CCradius+coaxialDLThickness,
			  CCradius+coaxialDLThickness};
    G4Polycone* coaxialDL =
      new G4Polycone("coaxialDL", 0.*deg, 360.*deg, 4, zSlice, rIn, rOut);

    G4LogicalVolume* coaxialDL_log =
      new G4LogicalVolume( coaxialDL, HpGe, "coaxialDL_log", 0, 0, 0 );
  
    // Place coaxial dead layer as a daugter of the leaf logical volume
    new G4PVPlacement( 0, G4ThreeVector(), coaxialDL_log, "coaxialDeadLayer",
		       leaf_log, false, 0);

    coaxialDL_log->SetVisAttributes(Vis_2);
  }
  
  // Cryosatat =================================================================

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

  cover_log->SetVisAttributes(Vis_2);

  // Copper backing ============================================================

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

  Cubox_log->SetVisAttributes(Vis_2);
 
  // Assemble the clover

  assemblyclover = new G4AssemblyVolume();

  Rot0.rotateZ(180.*deg);

  assemblyclover->AddPlacedVolume(leaf_log,Leaf0Pos,&Rot0);
 
  Rot0.rotateZ(90.*deg);

  assemblyclover->AddPlacedVolume(leaf_log,Leaf1Pos,&Rot0);

  Rot0.rotateZ(90.*deg);
  
  assemblyclover->AddPlacedVolume(leaf_log,Leaf2Pos,&Rot0);

  Rot0.rotateZ(90.*deg);

  assemblyclover->AddPlacedVolume(leaf_log,Leaf3Pos,&Rot0);

  assemblyclover->AddPlacedVolume(cover_log,coverpos,&wallrot);

  assemblyclover->AddPlacedVolume(Cubox_log,Cuboxpos,&wallrot);

  if(orientation == "left")
    DetCode = 31000 + 4*31 - 1;
  else if(orientation == "right")
    DetCode = 31000 + 4*32 - 1;

  assemblyclover->MakeImprint(expHall_log, DetPos, &DetRot, DetCode);

  return;
}
//---------------------------------------------------------------------
G4IntersectionSolid* Clover_Detector::Bulletized_Boxed_Cylinder(G4double offset)
{
  // Cylindrical crystal
  G4Tubs* detector = 
    new G4Tubs("detector", 0, Radius-offset, (Length-torusradius)/2.-offset,
	       0., 360.*deg);

  // Bulletized cap (height torusradius along z)
  G4Torus* torus = new G4Torus("torus", 0, torusradius - offset,
			       Radius - torusradius - offset,
			       0., 360.*deg);

  G4Box* torusbox = new G4Box("torusbox", Radius + torusradius*3., 
			      Radius + torusradius*3., 
			      torusradius - offset);
 
  G4Tubs* torustube =
    new G4Tubs("torustube", 0, Radius - torusradius - offset,
	       torusradius - offset, 
	       0., 360.*deg);

  G4UnionSolid* torus1 =
    new G4UnionSolid("torus2", torustube, torus,
		     G4Transform3D(G4RotationMatrix(),G4ThreeVector()));

  G4SubtractionSolid* torus2 =
    new G4SubtractionSolid("torus1", torus1, torusbox, 
			   G4Transform3D(G4RotationMatrix(),
					 G4ThreeVector(0.,0.,
						       torusradius-offset)));
  // Bulletized cylinder
  G4UnionSolid* bevel =
    new G4UnionSolid("bevel", detector, torus2,
		     G4Transform3D(G4RotationMatrix(),
				   G4ThreeVector(0,0,
						 -(Length-torusradius)/2.
						 +offset)));

  // ... with planar sides
  G4Box* box = new G4Box("box", boxlength/2.-offset, boxlength/2.-offset,
			 1.1*Length);

  G4IntersectionSolid* intersect =
    new G4IntersectionSolid ("intersect", bevel, box, 
			     G4Transform3D(G4RotationMatrix(),
					   G4ThreeVector(CCoffset, 
							 CCoffset, 
							 0.)));
  return intersect;
}
//---------------------------------------------------------------------
void Clover_Detector::MakeSensitive(TrackerGammaSD* TrackerGamma)
{
  leaf_log->SetSensitiveDetector(TrackerGamma);
}
//---------------------------------------------------------------------
