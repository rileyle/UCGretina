#ifdef LHTARGET
#include "Target_LH.hh"

Target::Target(Materials* mat)
{
  materials=mat;

  buildSled=false;

  targetCellType = "thick";

  lH2            = materials->FindMaterial("G4_lH2");
  vacuum         = materials->FindMaterial("G4_Galactic");
  Aluminum       = materials->FindMaterial("Al");
  Copper         = materials->FindMaterial("Cu");
  StainlessSteel = materials->FindMaterial("ssteel");
  Kapton         = materials->FindMaterial("Kapton");

  Pos = new G4ThreeVector(0.,0.,0.);
  NoRot = G4RotationMatrix::IDENTITY;
  TargetMaterial = lH2;
  targetDensity = 75.36*mg/cm3; // lH2 @ 16 K (NIST)
  targetAngle = 30.0*deg;       // NSCL GRETINA configuration

  // Beam tube
  BeamTubeRmax = 3.0*2.54*cm;
  BeamTubeRmin = BeamTubeRmax - 0.083*2.54*cm;
  BeamTubeDz   = 44.*cm; //LR (approx target to gate valve flange)

  BeamTeeRmax = 2.0*2.54*cm;
  BeamTeeRmin = BeamTeeRmax - 0.065*2.54*cm;
  BeamTeeDz   = 24.2*cm/2.;

  // Greta LH chamber (science fiction)
  GretaChamberRmin = 178.5*mm;
  GretaChamberRmax = 180.0*mm;
  GretaBTrmin = 18.5*mm;
  GretaBTrmax = 20.0*mm;
  GretaBTDz = (127.6/2*cm - GretaChamberRmax)/2.;

  BulgeDz = 1.0*mm;

  windowThickness = 125.*um;
  windows = false;

  NStep=20;
  
  sourceFrame = "";

  Greta = false;

  Cutaway = false;

}

Target::~Target()
{ delete target_limits;}

//-----------------------------------------------------------------------------
G4VPhysicalVolume* Target::Construct(G4LogicalVolume* experimentalHall_log)
{

  expHall_log=experimentalHall_log;

  LHTarget = new G4AssemblyVolume();

  if(Greta)
    ConstructGretaChamber();
  else
    ConstructBeamTube();

  // Target and Cell
  if(targetCellType == "thick")
    Target_phys = Construct200mgCell();
  else if(targetCellType == "thin")
    Target_phys = Construct50mgCell();
  else if(targetCellType == "empty")
    Target_phys = ConstructEmptyCell();
  else if(targetCellType == "notarget")
    Target_phys = ConstructNoCell();
  else
    G4cout << "Target_LH: unknown target cell type: "
	   << targetCellType << G4endl;

  if(targetCellType != "notarget")
    ConstructCryo();

  // Place the target system.
  G4RotationMatrix LHTargetRot = NoRot;
  LHTargetRot.rotateZ(targetAngle);

  LHTarget->MakeImprint(expHall_log, *Pos, &LHTargetRot);

  if(buildSled) BuildSled();

  return Target_phys;

}
//-----------------------------------------------------------------------------
void Target::ConstructBeamTube(){

  // Beam line segment

  G4Tubs* BeamTube0 = new G4Tubs("BeamTube0", BeamTubeRmin, BeamTubeRmax, 
				BeamTubeDz, 0.*deg, 360.*deg); 
  G4Tubs* BeamTee0 = new G4Tubs("BeamTee0", BeamTeeRmin, BeamTeeRmax, 
				BeamTeeDz, 0.*deg, 360.*deg); 

  G4Tubs* BeamTubeCutout = new G4Tubs("BeamTubeCutout", 0., BeamTeeRmax, 
				      BeamTeeDz, 0.*deg, 360.*deg); 

  G4RotationMatrix BeamTeeRot = NoRot;
  BeamTeeRot.rotateX(90.0*deg);

  G4Tubs* BeamTeeCutout = new G4Tubs("BeamTeeCutout", 0., BeamTubeRmin, 
				     BeamTubeDz, 0.*deg, 360.*deg); 

  G4SubtractionSolid* BeamTube 
    = new G4SubtractionSolid("BeamTube", BeamTube0, BeamTubeCutout,
			     G4Transform3D(BeamTeeRot,
					   G4ThreeVector(0., BeamTeeDz, 0.)));
  G4SubtractionSolid* BeamTee 
    = new G4SubtractionSolid ("BeamTee", BeamTee0, BeamTeeCutout,
			      G4Transform3D(BeamTeeRot,
					    G4ThreeVector(0., 0., BeamTeeDz)));

  G4LogicalVolume* BeamTube_log
    = new G4LogicalVolume(BeamTube, Aluminum, "BeamTube_log", 0, 0, 0);

  G4ThreeVector Origin = G4ThreeVector(0., 0., 0.);
  
  LHTarget->AddPlacedVolume(BeamTube_log, Origin, &NoRot);

  // Tee

  G4ThreeVector BeamTeeOffset = G4ThreeVector(0., BeamTeeDz, 0.);
  G4LogicalVolume* BeamTee_log
    = new G4LogicalVolume(BeamTee, Aluminum, "BeamTee_log", 0, 0, 0);

  LHTarget->AddPlacedVolume(BeamTee_log, BeamTeeOffset, &BeamTeeRot);

  G4Tubs* BeamTeeFlange = new G4Tubs("BeamTeeFlange", 
				     4.915*cm, 9.0*cm, 1.194*cm/2.,
				     0., 360.*deg);

  G4LogicalVolume* BeamTeeFlange_log
    = new G4LogicalVolume(BeamTeeFlange, Aluminum, "BeamTeeFlange_log", 0, 0, 0);

  G4RotationMatrix BeamTeeFlangeRot = NoRot;
  BeamTeeFlangeRot.rotateX(-90.0*deg);

  G4ThreeVector BeamTeeFlangePos = G4ThreeVector(0., 
						 24.2*cm + 1.194*cm/2.,
						 0.);

  LHTarget->AddPlacedVolume(BeamTeeFlange_log, 
			    BeamTeeFlangePos, 
			    &BeamTeeFlangeRot);

  // Bellows ===========================

  // DIMENSIONS APPROXIMATE
  G4Tubs* BellowsFlange = new G4Tubs("BellowsFlange",
				     4.*2.54/2*cm, 9.*2.54/2*cm, 
				     0.5*2.54/2*cm, 0.*deg, 360.*deg); 

  G4LogicalVolume* BellowsFlange_log
    = new G4LogicalVolume(BellowsFlange, StainlessSteel, "BellowsFlange_log", 0, 0, 0);

  G4RotationMatrix BellowsFlangeRot = NoRot;
  BellowsFlangeRot.rotateX(-90.0*deg);

  G4ThreeVector BellowsFlangePos = G4ThreeVector(0., 
						 24.2*cm + 1.194*cm
						 +0.75/2*cm,
						 0.);

  LHTarget->AddPlacedVolume(BellowsFlange_log, 
			    BellowsFlangePos, 
			    &BellowsFlangeRot);


  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0, 0.3); 
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(false);

  BeamTube_log->SetVisAttributes(Vis);
  BeamTee_log->SetVisAttributes(Vis);
  BeamTeeFlange_log->SetVisAttributes(Vis);
  BellowsFlange_log->SetVisAttributes(Vis);

}
//-----------------------------------------------------------------------------
void Target::ConstructGretaChamber(){

  G4double Phi0 = 0.*deg;
  G4double dPhi = 360.*deg;
  if(Cutaway){
    Phi0 = -90.*deg;
    dPhi = 180.*deg;
  }

  G4ThreeVector Origin = G4ThreeVector(0., 0., 0.);

  // Beam line segments

  G4Tubs* BeamTube = new G4Tubs("BeamTube", GretaBTrmin, GretaBTrmax, 
				GretaBTDz, Phi0, dPhi); 

  G4LogicalVolume* BeamTube_log = new G4LogicalVolume(BeamTube, 
						      Aluminum,
						      "BeamTube_log",
						      0, 0, 0);

  G4ThreeVector *GretaBTPos = new G4ThreeVector(0., 0., 0.);
  GretaBTPos->setZ(sqrt(GretaChamberRmin*GretaChamberRmin 
			- GretaBTrmin*GretaBTrmin) + GretaBTDz);

  LHTarget->AddPlacedVolume(BeamTube_log, *GretaBTPos, &NoRot);

  GretaBTPos->setZ(-sqrt(GretaChamberRmin*GretaChamberRmin 
			- GretaBTrmin*GretaBTrmin) - GretaBTDz);

  LHTarget->AddPlacedVolume(BeamTube_log, *GretaBTPos, &NoRot);

  // Target tee
  G4Tubs* BeamTee = new G4Tubs("BeamTee", BeamTeeRmin, BeamTeeRmax, 
			       BeamTeeDz, Phi0, dPhi); 

  G4RotationMatrix BeamTeeRot = NoRot;
  BeamTeeRot.rotateX(90.0*deg);

  G4ThreeVector* BeamTeeOffset 
    = new G4ThreeVector(0., sqrt(GretaChamberRmin*GretaChamberRmin 
				 -BeamTeeRmax*BeamTeeRmax) + BeamTeeDz, 
			0.);

  G4LogicalVolume* BeamTee_log
    = new G4LogicalVolume(BeamTee, Aluminum, "BeamTee_log", 0, 0, 0);

  LHTarget->AddPlacedVolume(BeamTee_log, *BeamTeeOffset, &BeamTeeRot);

  // Chamber
  G4Sphere* sphere = new G4Sphere("sphere", GretaChamberRmin, GretaChamberRmax, 
				  Phi0, dPhi, 0., 180.*deg);

  G4Tubs* drill = new G4Tubs("drill", 0, GretaBTrmax, 1.1*GretaChamberRmax, 
			     0., 360.*deg);

  G4SubtractionSolid* Chamber = new G4SubtractionSolid("Chamber", 
						       sphere, drill);

  G4Tubs* BeamTeeDrill = new G4Tubs("BeamTeeDrill", 0., BeamTeeRmax, 
				    BeamTeeDz, 0.*deg, 360.*deg); 

  Chamber = new G4SubtractionSolid("Chamber", 
				   Chamber, BeamTeeDrill,
				   G4Transform3D(BeamTeeRot,
						 G4ThreeVector(0., 
							       BeamTeeDz, 
							       0.)));

  G4LogicalVolume* Chamber_log 
    = new G4LogicalVolume(Chamber, Aluminum, "Chamber_log",
			  0, 0, 0);

  LHTarget->AddPlacedVolume(Chamber_log, Origin, &NoRot);

  // Flanges
  G4Tubs* BeamTeeFlange = new G4Tubs("BeamTeeFlange", 
				     4.915*cm, 9.0*cm, 1.194*cm/2.,
				     Phi0, dPhi);

  G4LogicalVolume* BeamTeeFlange_log
    = new G4LogicalVolume(BeamTeeFlange, Aluminum, "BeamTeeFlange_log", 0, 0, 0);

  G4RotationMatrix BeamTeeFlangeRot = NoRot;
  BeamTeeFlangeRot.rotateX(-90.0*deg);

  G4ThreeVector BeamTeeFlangePos 
    = G4ThreeVector(0., 
		    2.*BeamTeeDz +
		    sqrt(GretaChamberRmin*GretaChamberRmin 
			 -BeamTeeRmax*BeamTeeRmax));

  LHTarget->AddPlacedVolume(BeamTeeFlange_log, 
			    BeamTeeFlangePos, 
			    &BeamTeeFlangeRot);

  // DIMENSIONS APPROXIMATE
  G4Tubs* BellowsFlange = new G4Tubs("BellowsFlange",
				     4.*2.54/2*cm, 9.*2.54/2*cm, 
				     0.5*2.54/2*cm, Phi0, dPhi); 

  G4LogicalVolume* BellowsFlange_log
    = new G4LogicalVolume(BellowsFlange, StainlessSteel, "BellowsFlange_log", 0, 0, 0);

  G4RotationMatrix BellowsFlangeRot = NoRot;
  BellowsFlangeRot.rotateX(-90.0*deg);

  G4ThreeVector BellowsFlangePos 
    = G4ThreeVector(0., 
		    2.*BeamTeeDz +
		    sqrt(GretaChamberRmin*GretaChamberRmin 
			 -BeamTeeRmax*BeamTeeRmax)
		    +0.75/2*cm,
		    0.);

  LHTarget->AddPlacedVolume(BellowsFlange_log, 
			    BellowsFlangePos, 
			    &BellowsFlangeRot);

  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0, 0.3); 
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(false);

  BeamTube_log->SetVisAttributes(Vis);
  BeamTee_log->SetVisAttributes(Vis);
  Chamber_log->SetVisAttributes(Vis);
  BeamTeeFlange_log->SetVisAttributes(Vis);
  BellowsFlange_log->SetVisAttributes(Vis);

}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* Target::Construct200mgCell(){

  // Liquid Hydrogen
  //
  // (Thick target, 200mg/cm^2 nominal) 
  //
  // The window frames have an inner radius 1.9 cm over
  // 1.1 cm at either end of the cell. The region between the frames
  // is 1.5 cm wide and has an the inner  radius of 2.0 cm. The total
  // length of the assembled cell excluding clamp rings = 3.0 cm. This
  // is the thickness of the liquid hydrogen without window bulge.

  TargetDz = 3.0/2.*cm;
  BulgeR  = 1.9*cm;
  Target_thickness = 2*TargetDz + 2*BulgeDz;

  G4ThreeVector Pos0 = *Pos;
  
  const G4double TzPlane[6] = {0.,     0.75*cm, 0.75*cm, 2.25*cm, 
			       2.25*cm, 3.0*cm};
  const G4double TrInner[6] = {0.,     0.,      0.,      0,       
			       0,       0.};
  const G4double TrOuter[6] = {1.9*cm, 1.9*cm,  2.0*cm,  2.0*cm,  
			       1.9*cm,  1.9*cm};

  G4Polycone* Target0 = new G4Polycone("Target0",  0., 360.*deg, 
				       6, TzPlane, TrInner, TrOuter);
  G4Ellipsoid* Bulge = new G4Ellipsoid("Bulge",BulgeR,BulgeR,BulgeDz);
  G4UnionSolid* Target1 =
    new G4UnionSolid("target1", Target0, Bulge,
		     G4Transform3D(NoRot,
				   G4ThreeVector(0., 0., 0.)));

  aTarget = new G4UnionSolid("target1", Target1, Bulge,
			     G4Transform3D(NoRot,
					   G4ThreeVector(0., 0., 2*TargetDz)));

  // The standard lH2 density doesn't apply here, so we have to set it.
  G4String name=TargetMaterial->GetName();
  G4double Z=TargetMaterial->GetZ();
  G4double A=TargetMaterial->GetA();
  TargetMaterial=new G4Material(name, Z,A,targetDensity);

  Target_log = new G4LogicalVolume(aTarget,TargetMaterial,"Target_log",0,0,0);
  target_limits= new G4UserLimits();
  target_limits->SetMaxAllowedStep(Target_thickness/NStep);
  Target_log->SetUserLimits(target_limits);

  G4Colour blue (0.0,0.0,1.0); 
  G4VisAttributes* Vis = new G4VisAttributes(blue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(true);

  Target_log->SetVisAttributes(Vis);

  G4ThreeVector TargetPos = G4ThreeVector(0., 0., -TargetDz);
  TargetPos += Pos0;

  Target_phys = new G4PVPlacement(G4Transform3D(NoRot,TargetPos),
				  Target_log,"Target",expHall_log,false,0);

  // Cell windows
  //
  // Include "Target" in the volume name so that it is included
  // in determining when the reaction product leaves the target.
  if(windows){
    G4SubtractionSolid* UpstreamWindow 
      = new G4SubtractionSolid("UpStrTargetWindow", Bulge, Bulge,
			       G4Transform3D(NoRot,
					     G4ThreeVector(0., 
							   0., 
							   windowThickness)));

    G4LogicalVolume* UpstreamWindow_log
      = new G4LogicalVolume(UpstreamWindow, Kapton, "UpStrTargetWindow_log", 
			    0, 0, 0);
    G4ThreeVector USwindowPos = G4ThreeVector(0., 0.,
					      -TargetDz-windowThickness);
    
    LHTarget->AddPlacedVolume(UpstreamWindow_log, USwindowPos, &NoRot);

    G4SubtractionSolid* DnstreamWindow 
      = new G4SubtractionSolid("DnStrTargetWindow", Bulge, Bulge,
			       G4Transform3D(NoRot,
					     G4ThreeVector(0., 
							   0., 
							   -windowThickness)));
    G4LogicalVolume* DnstreamWindow_log
      = new G4LogicalVolume(DnstreamWindow, Kapton, "DnStrTargetWindow_log", 
			    0, 0, 0);
    G4ThreeVector DSwindowPos = G4ThreeVector(0., 0., 
					      TargetDz+windowThickness);
    
    LHTarget->AddPlacedVolume(DnstreamWindow_log, DSwindowPos, &NoRot);

  }
#if(1)
  // Target Cell
  //
  // The window frames and clamp rings have inner radius 1.9 cm over
  // 1.1 cm at either end of the cell. The region between the frames
  // is 1.5 cm wide and has an the inner  radius of 2.0 cm. The total
  // length of the assembled cell including clamp rings is 3.7 cm.

  const G4double CzPlane[6] = {0.,     1.1*cm, 1.1*cm, 2.6*cm, 2.6*cm, 3.7*cm};
  const G4double CrInner[6] = {1.9*cm, 1.9*cm, 2.0*cm, 2.0*cm, 1.9*cm, 1.9*cm};
  const G4double CrOuter[6] = {2.7*cm, 2.7*cm, 2.7*cm, 2.7*cm, 2.7*cm, 2.7*cm};

  G4Polycone* TargetCell = new G4Polycone("TargetCell",  0., 360.*deg, 
					   6, CzPlane, CrInner, CrOuter);

  G4LogicalVolume* TargetCell_log
    = new G4LogicalVolume(TargetCell, Aluminum, "TargetCell_log", 0, 0, 0);

  G4ThreeVector TargetCellPos = G4ThreeVector(0., 0., -3.7*cm/2.);

  LHTarget->AddPlacedVolume(TargetCell_log, TargetCellPos, &NoRot);
#endif
  return Target_phys;

}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* Target::Construct50mgCell(){

  // Liquid Hydrogen
  //
  // (Thick target, 50 mg/cm^2 nominal) 
  //
  // The window frames have an inner radius 1.59 cm over
  // 1.15 cm at either end of the cell. The region between the frames
  // is 0.695 cm wide and has an inner radius of 2.0 cm. This is
  // the thickness of the liquid hydrogen without the bulges.
  // The total length of the assembled cell is 3.0 cm.

  TargetDz = 0.695*cm/2.;
  BulgeR  = 1.59*cm;
  Target_thickness = 2*TargetDz + 2*BulgeDz;

  G4Tubs* Target0 = new G4Tubs("Target0", 0., 2.0*cm, TargetDz,
			       0., 360.*deg);
			       
  G4Ellipsoid* Bulge = new G4Ellipsoid("Bulge",BulgeR,BulgeR,BulgeDz);
  G4UnionSolid* Target1 =
    new G4UnionSolid("target1", Target0, Bulge,
		     G4Transform3D(NoRot,
				   G4ThreeVector(0., 0., -TargetDz)));

  aTarget = new G4UnionSolid("target1", Target1, Bulge,
			     G4Transform3D(NoRot,
					   G4ThreeVector(0., 0., TargetDz)));

  Target_log = new G4LogicalVolume(aTarget,lH2,"Target_log",0,0,0);
  target_limits= new G4UserLimits();
  target_limits->SetMaxAllowedStep(Target_thickness/NStep);
  Target_log->SetUserLimits(target_limits);

  G4Colour blue (0.0,0.0,1.0); 
  G4VisAttributes* Vis = new G4VisAttributes(blue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(true);

  Target_log->SetVisAttributes(Vis);

  Target_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),
				  Target_log,"Target",expHall_log,false,0);

  // Cell windows
  if(windows){
    G4SubtractionSolid* UpstreamWindow 
      = new G4SubtractionSolid("UpstreamWindow", Bulge, Bulge,
			       G4Transform3D(NoRot,
					     G4ThreeVector(0., 
							   0., 
							   windowThickness)));

    G4LogicalVolume* UpstreamWindow_log
      = new G4LogicalVolume(UpstreamWindow, Kapton, "UpstreamWindow_log", 
			    0, 0, 0);

    G4Colour red (1.0,0.0,0.0); 
    G4VisAttributes* Vis2 = new G4VisAttributes(red);
    Vis2->SetVisibility(true);
    UpstreamWindow_log->SetVisAttributes(Vis2);

    G4ThreeVector USwindowPos = G4ThreeVector(0., 0.,
					      -TargetDz-windowThickness);
    
    LHTarget->AddPlacedVolume(UpstreamWindow_log, USwindowPos, &NoRot);

    G4SubtractionSolid* DnstreamWindow 
      = new G4SubtractionSolid("DnstreamWindow", Bulge, Bulge,
			       G4Transform3D(NoRot,
					     G4ThreeVector(0., 
							   0., 
							   -windowThickness)));
    G4LogicalVolume* DnstreamWindow_log
      = new G4LogicalVolume(DnstreamWindow, Kapton, "DnstreamWindow_log", 
			    0, 0, 0);
    DnstreamWindow_log->SetVisAttributes(Vis2);
    G4ThreeVector DSwindowPos = G4ThreeVector(0., 0., 
					      TargetDz+windowThickness);

    LHTarget->AddPlacedVolume(DnstreamWindow_log, DSwindowPos, &NoRot);

  }
#if(1)
  // Target Cell
  //
  // The window frames and clamp rings have inner radius 1.59 cm over
  // 1.15 cm at either end of the cell. The region between the frames
  // is 0.695 cm wide and has an inner radius of 2.0 cm. The total
  // length of the assembled cell is 3.0 cm.

  const G4double CzPlane[6] = {0.,      1.15*cm, 1.15*cm, 1.85*cm, 
			       1.85*cm, 3.0*cm};
  const G4double CrInner[6] = {1.59*cm, 1.59*cm, 2.0*cm,  2.0*cm,
			       1.59*cm, 1.59*cm};
  const G4double CrOuter[6] = {2.7*cm, 2.7*cm,  2.7*cm,  2.7*cm,  
			       2.7*cm,  2.7*cm};


  G4Polycone* TargetCell = new G4Polycone("TargetCell",  0., 360.*deg, 
					   6, CzPlane, CrInner, CrOuter);

  G4LogicalVolume* TargetCell_log
    = new G4LogicalVolume(TargetCell, Aluminum, "TargetCell_log", 0, 0, 0);

  G4ThreeVector TargetCellPos = G4ThreeVector(0., 0., -3.0*cm/2.);
  
  LHTarget->AddPlacedVolume(TargetCell_log, TargetCellPos, &NoRot);
#endif
  return Target_phys;

}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* Target::ConstructEmptyCell(){

  // Set these to dummy values. They aren't used in this context.
  TargetDz = 1.*cm/2.;
  BulgeR  = 2.*cm;
  Target_thickness = TargetDz + 2*BulgeDz;

  G4ThreeVector Origin = G4ThreeVector(0., 0., 0.);
  
  // G4_Galactic target volume (the code needs a G4PhysicalVolume,
  // even for source simulations).

  G4Tubs* Target0 = new G4Tubs("Target0", 0., 2.0*cm, TargetDz,
			       0., 360.*deg);
			       
  G4Ellipsoid* Bulge = new G4Ellipsoid("Bulge",BulgeR,BulgeR,BulgeDz);
  G4UnionSolid* Target1 =
    new G4UnionSolid("target1", Target0, Bulge,
		     G4Transform3D(NoRot,
				   G4ThreeVector(0., 0., -TargetDz)));

  aTarget = new G4UnionSolid("target1", Target1, Bulge,
			     G4Transform3D(NoRot,
					   G4ThreeVector(0., 0., TargetDz)));

  Target_log = new G4LogicalVolume(aTarget,vacuum,"Target_log",0,0,0);
  target_limits= new G4UserLimits();
  target_limits->SetMaxAllowedStep(Target_thickness/NStep);
  Target_log->SetUserLimits(target_limits);

  G4Colour blue (0.0,0.0,1.0); 
  G4VisAttributes* Vis = new G4VisAttributes(blue);
  Vis->SetVisibility(false);

  Target_log->SetVisAttributes(Vis);

  Target_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),
				  Target_log,"Target",expHall_log,false,0);

  // Target Cell Body

  G4Tubs* TargetCell = new G4Tubs("TargetCell", 2.0*cm, 2.7*cm, 2.3*cm/2.,
				  0., 360.*deg);

  G4LogicalVolume* TargetCell_log
    = new G4LogicalVolume(TargetCell, Aluminum, "TargetCell_log", 0, 0, 0);

  //  LHTarget->AddPlacedVolume(TargetCell_log, *Pos, &NoRot);
  LHTarget->AddPlacedVolume(TargetCell_log, Origin, &NoRot);

  return Target_phys;

}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* Target::ConstructNoCell(){

  // Set these to dummy values. They aren't used in this context.
  TargetDz = 1.*cm/2.;
  BulgeR  = 2.*cm;
  Target_thickness = TargetDz + 2*BulgeDz;

  // G4_Galactic target volume (the code needs a G4PhysicalVolume,
  // even for source simulations).

  G4Tubs* Target0 = new G4Tubs("Target0", 0., 2.0*cm, TargetDz,
			       0., 360.*deg);
			       
  G4Ellipsoid* Bulge = new G4Ellipsoid("Bulge",BulgeR,BulgeR,BulgeDz);
  G4UnionSolid* Target1 =
    new G4UnionSolid("target1", Target0, Bulge,
		     G4Transform3D(NoRot,
				   G4ThreeVector(0., 0., -TargetDz)));

  aTarget = new G4UnionSolid("target1", Target1, Bulge,
			     G4Transform3D(NoRot,
					   G4ThreeVector(0., 0., TargetDz)));

  Target_log = new G4LogicalVolume(aTarget,vacuum,"Target_log",0,0,0);
  target_limits= new G4UserLimits();
  target_limits->SetMaxAllowedStep(Target_thickness/NStep);
  Target_log->SetUserLimits(target_limits);

  G4Colour blue (0.0,0.0,1.0); 
  G4VisAttributes* Vis = new G4VisAttributes(blue);
  Vis->SetVisibility(false);

  Target_log->SetVisAttributes(Vis);

  Target_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),
				  Target_log,"Target",expHall_log,false,0);

  return Target_phys;

}
//-----------------------------------------------------------------------------
void Target::ConstructCryo(){
  // Target Cell Stem and Flange =======

  G4ThreeVector Pos0 = *Pos;
  
  const G4double CSzPlane[4] = {0.,      1.69*cm, 1.69*cm, 2.32*cm};
  const G4double CSrInner[4] = {0.16*cm, 0.16*cm, 0.16*cm, 0.16*cm};
  const G4double CSrOuter[4] = {0.48*cm, 0.48*cm, 1.43*cm, 1.43*cm};

  G4Polycone* TargetStem0 = new G4Polycone("TargetStem0",  0., 360.*deg, 
					   4, CSzPlane, CSrInner, CSrOuter);

  G4RotationMatrix TargetStemRot = NoRot;
  TargetStemRot.rotateX(-90.0*deg);
  G4Tubs* TargetStemCutout = new G4Tubs("TargetStemCutout", 0., 2.7*cm, 1*cm, 
					0., 360.*deg);
  G4SubtractionSolid* TargetStem
    = new G4SubtractionSolid("TargetStem", TargetStem0, TargetStemCutout,
			     G4Transform3D(TargetStemRot,
					   G4ThreeVector(0., 0., -2.7*cm)));

  G4LogicalVolume* TargetStem_log
    = new G4LogicalVolume(TargetStem, Aluminum, "TargetStem_log", 0, 0, 0);

  G4double CellandStemHeight = 4.98*cm;

  G4ThreeVector TargetStemPos = G4ThreeVector(0., 2.65*cm, 0.);

  LHTarget->AddPlacedVolume(TargetStem_log, TargetStemPos, &TargetStemRot);

  // Cryocooler Extension Bar ==========
  G4double ExtensionBarHeight = 19.05*cm;

  const G4double EBzPlane[9] = {0.,       0.63*cm,  0.63*cm,  2.63*cm, 4.44*cm, 
				4.44*cm, 18.41*cm, 18.41*cm, 19.50*cm};
  const G4double EBrInner[9] = {0.48*cm, 0.48*cm,   0.48*cm,  0.48*cm,  0.,
				0.,      0.,        0.,       0.};
  const G4double EBrOuter[9] = {1.42*cm, 1.42*cm,  0.71*cm,   0.71*cm, 0.71*cm, 
				1.75*cm, 1.75*cm,  2.60*cm,   2.60*cm};

  G4Polycone* ExtensionBar = new G4Polycone("ExtensionBar",  0., 360.*deg, 
					    9, EBzPlane, EBrInner, EBrOuter);

  G4RotationMatrix ExtensionBarRot = NoRot;
  ExtensionBarRot.rotateX(-90.0*deg);

  G4LogicalVolume* ExtensionBar_log
    = new G4LogicalVolume(ExtensionBar, Copper, "ExtensionBar_log", 0, 0, 0);

  G4ThreeVector ExtensionBarPos = G4ThreeVector(0., CellandStemHeight, 0.);
  
  LHTarget->AddPlacedVolume(ExtensionBar_log, ExtensionBarPos, 
			    &ExtensionBarRot);

  // Heater Block ======================
  G4double HeaterBlockHeight = 1.27*cm;
  G4Tubs* HeaterBlock = new G4Tubs("HeaterBlock", 0., 
				   2.6*cm, HeaterBlockHeight/2., 
				   0., 360.*deg);

  G4LogicalVolume* HeaterBlock_log
    = new G4LogicalVolume(HeaterBlock, Copper, "HeaterBlock_log", 0, 0, 0);

  G4ThreeVector HeaterBlockPos = G4ThreeVector(0., 
					       CellandStemHeight + 
					       ExtensionBarHeight +
					       HeaterBlockHeight/2.,
					       0.);
  
  LHTarget->AddPlacedVolume(HeaterBlock_log, HeaterBlockPos, &TargetStemRot);

  // Cryocooler 2nd Stage ==============
  //  G4double Cryocooler2ndStageHeight = 20.68*cm;
  const G4double S2zPlane[10] = {0.,        0.61*cm,  0.61*cm,  6.13*cm,  
				 6.13*cm,  17.02*cm, 17.02*cm, 20.02*cm, 
				 20.02*cm, 20.68*cm};
  const G4double S2rInner[10] = {0.,        0.,       0.,       0.,
				 0.,        0.,       0.,       0.,
				 0.,        0.};
  const G4double S2rOuter[10] = {2.6*cm,    2.6*cm,   1.75*cm,  1.75*cm,  
				 1.39*cm,   1.39*cm,  3.*cm,    3.*cm,
				 4.51*cm,   4.51*cm};

  G4Polycone* Cryocooler2ndStage = new G4Polycone("Cryocooler2ndStage",  
						  0., 360.*deg,  
						  10, 
						  S2zPlane, 
						  S2rInner, 
						  S2rOuter);

  G4RotationMatrix Cryocooler2ndStageRot = NoRot;
  Cryocooler2ndStageRot.rotateX(-90.0*deg);

  G4LogicalVolume* Cryocooler2ndStage_log
    = new G4LogicalVolume(Cryocooler2ndStage, Copper, 
			  "Cryocooler2ndStage_log", 0, 0, 0);

  G4ThreeVector Cryocooler2ndStagePos = G4ThreeVector(0., 
						      CellandStemHeight +
						      ExtensionBarHeight +
						      HeaterBlockHeight,
						      0.);
  
  LHTarget->AddPlacedVolume(Cryocooler2ndStage_log, 
			    Cryocooler2ndStagePos, 
			    &Cryocooler2ndStageRot);

  // Gas Line ==========================

  G4double GasLineRmin  =  2.46*mm;
  G4double GasLineRmax  =  3.17*mm;
  G4double GasLineRbend = 14.30*mm;

  // Gas Line: 0.1 inch straight segment attached to the extension bar

  G4Tubs* GasLine1 = new G4Tubs("GasLine1", GasLineRmin, GasLineRmax,
				2.54*mm/2., 0., 360.*deg);

  G4LogicalVolume* GasLine1_log
    = new G4LogicalVolume(GasLine1, StainlessSteel, "GasLine1_log", 0, 0, 0);

  G4ThreeVector GasLine1Pos = G4ThreeVector(0., 75.2*mm, 7.11*mm + 2.54*mm/2.);
  
  LHTarget->AddPlacedVolume(GasLine1_log, GasLine1Pos, &NoRot);

  // Gas Line: First bend (180 deg)

  G4Torus* GasLine2 = new G4Torus("GasLine2", GasLineRmin, GasLineRmax,
				  GasLineRbend, 180.*deg, 180*deg);

  G4LogicalVolume* GasLine2_log
    = new G4LogicalVolume(GasLine2, StainlessSteel, "GasLine2_log", 0, 0, 0);

  G4RotationMatrix GasLine2Rot = NoRot;
  GasLine2Rot.rotateX(-90.0*deg);

  G4ThreeVector GasLine2Pos = G4ThreeVector(GasLineRbend, 
					    75.2*mm, 
					    7.11*mm + 2.54*mm/2.);
  
  LHTarget->AddPlacedVolume(GasLine2_log, GasLine2Pos, &GasLine2Rot);

  // Gas Line: 0.75 inch straight segment

  G4Tubs* GasLine3 = new G4Tubs("GasLine3", GasLineRmin, GasLineRmax,
				19.05*mm/2., 0., 360.*deg);

  G4LogicalVolume* GasLine3_log
    = new G4LogicalVolume(GasLine3, StainlessSteel, "GasLine3_log", 0, 0, 0);

  G4ThreeVector GasLine3Pos = G4ThreeVector(2*GasLineRbend, 75.2*mm, 0.);
  
  LHTarget->AddPlacedVolume(GasLine3_log, GasLine3Pos, &NoRot);

  // Gas Line: Second bend (121 deg)

  G4Torus* GasLine4 = new G4Torus("GasLine4", GasLineRmin, GasLineRmax,
				  GasLineRbend, 0.*deg, 121.*deg);

  G4LogicalVolume* GasLine4_log
    = new G4LogicalVolume(GasLine4, StainlessSteel, "GasLine4_log", 0, 0, 0);

  G4RotationMatrix GasLine4Rot = NoRot;
  GasLine4Rot.rotateX(-90.0*deg);

  G4ThreeVector GasLine4Pos = G4ThreeVector(GasLineRbend, 
					    75.2*mm, 
					    -7.11*mm - 2.54*mm/2.);
  
  LHTarget->AddPlacedVolume(GasLine4_log, GasLine4Pos, &GasLine2Rot);

  // Gas Line: 1.098 inch straight segment

  G4Tubs* GasLine5 = new G4Tubs("GasLine5", GasLineRmin, GasLineRmax,
				27.89*mm/2., 0., 360.*deg);

  G4LogicalVolume* GasLine5_log
    = new G4LogicalVolume(GasLine5, StainlessSteel, "GasLine5_log", 0, 0, 0);

  G4RotationMatrix GasLine5Rot = NoRot;
  GasLine5Rot.rotateY(-59.0*deg);

  G4ThreeVector GasLine5Pos = G4ThreeVector(-5.1, 75.2*mm, -13.45);
  
  LHTarget->AddPlacedVolume(GasLine5_log, GasLine5Pos, &GasLine5Rot);
 
  // Third bend (90 deg to vertical)

  G4Torus* GasLine6 = new G4Torus("GasLine6", GasLineRmin, GasLineRmax,
				  GasLineRbend, 0.*deg, 90.*deg);

  G4LogicalVolume* GasLine6_log
    = new G4LogicalVolume(GasLine6, StainlessSteel, "GasLine6_log", 0, 0, 0);

  G4RotationMatrix GasLine6Rot = NoRot;
  GasLine6Rot.rotateZ(-90.0*deg);
  GasLine6Rot.rotateY(-149.0*deg);

  G4ThreeVector GasLine6Pos = G4ThreeVector(-17.1*mm, 
					    75.2*mm+GasLineRbend, 
					    -6.2*mm);
  
  LHTarget->AddPlacedVolume(GasLine6_log, GasLine6Pos, &GasLine6Rot);

  // Gas Line: 8.25 inch vertical straight segment

  G4Tubs* GasLine7 = new G4Tubs("GasLine7", GasLineRmin, GasLineRmax,
				209.55*mm/2., 0., 360.*deg);

  G4LogicalVolume* GasLine7_log
    = new G4LogicalVolume(GasLine7, StainlessSteel, "GasLine7_log", 0, 0, 0);

  G4RotationMatrix GasLine7Rot = NoRot;
  GasLine7Rot.rotateX(90.0*deg);

  G4ThreeVector GasLine7Pos = G4ThreeVector(-GasLineRbend - 15.*mm, 
					    75.2*mm + 209.55*mm/2.
					    + GasLineRbend, 
					    1.2);
  
  LHTarget->AddPlacedVolume(GasLine7_log, GasLine7Pos, &GasLine7Rot);

  // Radiation Shield ==================

  G4double Phi0 = 0.*deg;
  G4double dPhi = 360.*deg;
  if(Cutaway){
    Phi0 = -90.*deg;
    dPhi = 180.*deg;
  }

  const G4double RSzPlane[6] = { 0.,       0.32*cm,  0.32*cm,  
				 49.19*cm, 49.19*cm, 50.14*cm};
  const G4double RSrInner[6] = { 0.,       0.,       3.64*cm,
				 3.64*cm,  3.64*cm,  3.64*cm};
  const G4double RSrOuter[6] = { 3.81*cm,  3.81*cm,  3.81*cm,
				 3.81*cm,  5.87*cm,  5.87*cm};

  G4Polycone* RadiationShield0 = new G4Polycone("RadiationShield0",  
						Phi0, dPhi,
						6, 
						RSzPlane, 
						RSrInner, 
						RSrOuter);

  G4RotationMatrix RadiationShieldRot = NoRot;
  RadiationShieldRot.rotateX(-90.0*deg);

  G4ThreeVector RadiationShieldPos = G4ThreeVector(0., 
						   -4.85*cm,
						   0.);
  
  G4Tubs* RadiationShieldCutout = new G4Tubs("RadiationShieldCutout", 
					     0., 2.22*cm, 5.0*cm, 
					     0., 360.*deg);

  G4SubtractionSolid* RadiationShield
    = new G4SubtractionSolid ("RadiationShield", 
			      RadiationShield0, RadiationShieldCutout,
			      G4Transform3D(RadiationShieldRot,
					    G4ThreeVector(0., 0., 4.85*cm)));

  G4LogicalVolume* RadiationShield_log
    = new G4LogicalVolume(RadiationShield, Aluminum, 
			  "RadiationShield_log", 0, 0, 0);

  LHTarget->AddPlacedVolume(RadiationShield_log, 
			    RadiationShieldPos, 
			    &RadiationShieldRot);

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
  G4cout<<"----> Using the " << targetCellType << " target cell." << G4endl;
  if(Greta)
      G4cout<<"----> Constructing the GRETA target." << G4endl;
  if(Cutaway)
      G4cout<<"----> Constructing the cutaway view. For visualization only!" << G4endl;
  G4cout<<"----> Target bulge thickness set to " << BulgeDz/mm << " mm." << G4endl;
  G4cout<<"----> Target angle set to " << targetAngle/deg << " deg." << G4endl;
  G4cout<<"----> Target material set to  "<<Target_log->GetMaterial()->GetName()<< G4endl;   
  G4cout<<"----> Target density:         "<<Target_log->GetMaterial()->GetDensity()/mg*cm3 << " mg/cm^3" << G4endl;   
  G4cout<<"----> Number of simulation steps in the target is set to "<<NStep<<G4endl;
}
//---------------------------------------------------------------------
void Target::setMaterial(G4String materialName)
{
  if(materialName == "G4_Galactic" ||
     materialName == "vacuum"){

    Target_log->SetMaterial(vacuum);
    G4cout<<"----> Target material set to     "<<Target_log->GetMaterial()->GetName()<< G4endl;

    G4Colour blue (0.0,0.0,1.0); 
    G4VisAttributes* Vis = new G4VisAttributes(blue);
    Vis->SetVisibility(false);

    Target_log->SetVisAttributes(Vis);

  } else {
    G4cout<<"----> Warning: User can't set the material of the LH Target to anything other than G4_Galactic. Ignoring Target::setMaterial("<< materialName << ")."<<G4endl;
  }
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
  Pos->setZ(d);
  
  G4cout <<"----> Target Z position is set to "
	 << G4BestUnit(Pos->getZ(), "Length")
	 << G4endl;
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
void Target::BuildSled()
{
  if(targetCellType != "notarget")
    G4cout<<"----> Warning: target sled specified with LH target. Proceeding with no sled. "<< G4endl;
  else{
    sledMaterial = materials->FindMaterial("G10");

    G4double tolerance = 0.1*mm;

    sledFrameThickness  = 0.02*2.54*cm;
    sledFrameRmin       = 2.5*2.54/2.0*cm;
    sledFrameRmax       = 3.0*2.54/2.0*cm;

    sledRunnerThickness = 0.031*2.54*cm;
    sledRunnerLength    = 4.5*2.54*cm;
    sledRunnerHeight    = 1.5*2.54*cm;
    sledRunnerRmin      = 1.3*2.54*cm;
    sledRunnerRmax      = 2.942*2.54*cm;

    sledBarThickness = 0.25*2.54*cm;
    sledBarLength    = 3.375*2.54*cm;
    sledBarLA        = 0.5*2.54*cm;
    sledBarHA        = 0.75*2.54*cm;
    sledBarLB        = sledBarLength - 2*sledBarLA;
    sledBarHBC        = 0.25*2.54*cm;
    sledBarLC        = 0.25*2.54*cm;

    sledFrame = new G4Tubs("sledFrame", sledFrameRmin, sledFrameRmax, sledFrameThickness/2.0, 0., 360.*deg);

    sledFrame_log = new G4LogicalVolume(sledFrame, sledMaterial, "sledFrame_log", 0, 0, 0);

    G4ThreeVector   *Pos0 = new G4ThreeVector(0., 0., 
					      Target_thickness/2.0 
					      + sledFrameThickness/2.0 
					      + tolerance);
    sledFrame_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), sledFrame_log, "sledFrame_phys", expHall_log, false, 0);

    sledRunnerBox  = new G4Box("sledRunnerBox", sledRunnerLength/2.0, sledRunnerHeight/2.0, sledRunnerThickness/2.0);
    sledRunnerTubs = new G4Tubs("sledRunnerTubs", sledRunnerRmin, sledRunnerRmax, sledRunnerThickness, 0., 360.*deg);
    Pos0->set(0., 1.5*2.54*cm, 0.);
    sledRunner     = new G4IntersectionSolid("sledRunner", sledRunnerBox, sledRunnerTubs, G4Transform3D(NoRot,*Pos0));

    sledRunner_log = new G4LogicalVolume(sledRunner, sledMaterial, "sledRunner_log", 0, 0, 0);

    Pos0->set(0., -1.5*2.54*cm, 
	      -sledBarLength/2.0 
	      - sledRunnerThickness/2.0 
	      + Target_thickness/2.0 
	      + sledFrameThickness/2.0);
    sledRunner1_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), sledRunner_log, "sledRunner1_phys", expHall_log, false, 0);

    Pos0->set(0., -1.5*2.54*cm, 
	      sledBarLength/2.0 
	      + sledRunnerThickness/2.0 
	      + Target_thickness/2.0 
	      + sledFrameThickness/2.0
	      + 2.*tolerance);
    sledRunner2_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), sledRunner_log, "sledRunner2_phys", expHall_log, false, 0);

    sledBarBoxA    = new G4Box("sledBarBoxA", sledBarThickness/2.0, sledBarHA/2.0,  sledBarLA/2.0);
    sledBarBoxB    = new G4Box("sledBarBoxB", sledBarThickness/2.0, sledBarHBC/2.0, sledBarLB/2.0);
    sledBarBoxC    = new G4Box("sledBarBoxC", sledBarThickness/2.0, sledBarHBC/2.0, sledBarLC/2.0);

    sledBarBoxA_log = new G4LogicalVolume(sledBarBoxA, sledMaterial, "sledBarBoxA_log", 0, 0, 0);
    sledBarBoxB_log = new G4LogicalVolume(sledBarBoxB, sledMaterial, "sledBarBoxB_log", 0, 0, 0);
    sledBarBoxC_log = new G4LogicalVolume(sledBarBoxC, sledMaterial, "sledBarBoxC_log", 0, 0, 0);

    sledBar        = new G4AssemblyVolume();
  
    sledBarZ1        = -sledBarLength/2.0 + 0.25*2.54*cm;
    sledBarY1        = 0.;
    Pos0 = new G4ThreeVector(0., sledBarY1, sledBarZ1-tolerance);
    sledBar->AddPlacedVolume(sledBarBoxA_log, *Pos0, &NoRot);

    sledBarZ2        = 0.;
    sledBarY2        = 0.;
    Pos0->set(0., sledBarY2, sledBarZ2);
    sledBar->AddPlacedVolume(sledBarBoxB_log, *Pos0, &NoRot);

    sledBarZ3        = -sledFrameThickness/2.0 - sledBarLC/2.0;
    sledBarY3        = 0.25*2.54*cm;
    Pos0->set(0., sledBarY3+tolerance, sledBarZ3-tolerance);
    sledBar->AddPlacedVolume(sledBarBoxC_log, *Pos0, &NoRot);

    sledBarZ4        = sledFrameThickness/2.0 + sledBarLC/2.0;
    sledBarY4        = sledBarY3;
    Pos0->set(0., sledBarY4+tolerance, sledBarZ4+tolerance);
    sledBar->AddPlacedVolume(sledBarBoxC_log, *Pos0, &NoRot);

    sledBarZ5        = sledBarLength/2.0 - 0.25*2.54*cm;
    sledBarY5        = 0.;
    Pos0->set(0., sledBarY5, sledBarZ5+tolerance);
    sledBar->AddPlacedVolume(sledBarBoxA_log, *Pos0, &NoRot);

    Pos0->set(-(1.011 + 1.287)/2.0*2.54*cm, -(1.011 + 1.287)/2.0*2.54*cm,
	      Target_thickness/2.0 
	      + sledFrameThickness/2.0 
	      + tolerance);
    *Pos0 += *Pos;
    G4RotationMatrix Rot0 = G4RotationMatrix::IDENTITY;  
    Rot0.rotateZ(-45.*deg);
    sledBar->MakeImprint(expHall_log, *Pos0, &Rot0);

    Pos0->set( (1.011 + 1.287)/2.0*2.54*cm, -(1.011 + 1.287)/2.0*2.54*cm,
	       Target_thickness/2.0 
	       + sledFrameThickness/2.0
	       + tolerance);
    *Pos0 += *Pos;
    Rot0 = G4RotationMatrix::IDENTITY;  
    Rot0.rotateZ(45.*deg);
    sledBar->MakeImprint(expHall_log, *Pos0, &Rot0);

    G4cout << "----> Including the target sled." << G4endl;
  }

}
//-------------------------------------------------------------------
#endif
