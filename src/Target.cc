#ifndef LHTARGET
#include "Target.hh"

Target::Target()
{
  buildSled=false;

  Target_side_x=50*mm;
  Target_side_y=50*mm;
  Target_thickness=0.1*mm;
  TargetMaterialName = "G4_Galactic";
  Target_density_scale = 1.0;
  Pos = new G4ThreeVector(0.,0.,0.);
  Rot = G4RotationMatrix::IDENTITY;
  Rot.rotateZ(45.*deg);
  NoRot = G4RotationMatrix::IDENTITY;
  NStep=20;

  sourceFrame = "";
}

Target::~Target()
{ delete target_limits;}

//-----------------------------------------------------------------------------
G4VPhysicalVolume* Target::Construct(G4LogicalVolume* experimentalHall_log)
{
  expHall_log=experimentalHall_log;

  TargetMaterial = G4Material::GetMaterial(TargetMaterialName);
  G4String name=TargetMaterial->GetName();
  G4double Z=TargetMaterial->GetZ();
  G4double A=TargetMaterial->GetA();
  G4double density=TargetMaterial->GetDensity();
  density*=Target_density_scale;
  TargetMaterial=new G4Material(name, Z,A,density);

  aTarget = new G4Box("target",Target_side_x/2.,Target_side_y/2.,Target_thickness/2.);

  Target_log = new G4LogicalVolume(aTarget,TargetMaterial,"Target_log",0,0,0);
  target_limits= new G4UserLimits();
  target_limits->SetMaxAllowedStep(Target_thickness/NStep);
  Target_log->SetUserLimits(target_limits);

  Target_phys = new G4PVPlacement(G4Transform3D(Rot,*Pos),Target_log,"Target",expHall_log,false,0);

  G4Colour blue (0.0, 0.0, 1.0); 
  G4VisAttributes* Vis_6 = new G4VisAttributes(blue);
  if(TargetMaterial->GetName() == "G4_Galactic"){
    Vis_6->SetVisibility(false);
  } else {
    Vis_6->SetVisibility(true);
    Vis_6->SetForceSolid(true);
  }

  Target_log->SetVisAttributes(Vis_6);

  if(buildSled) BuildSled();

  return Target_phys;
}

//-----------------------------------------------------------------------------
void Target::setX(G4double X)
{

  Target_side_x=X;
  G4cout<<"----> Target side X is set to "<<G4BestUnit(Target_side_x,"Length")<<G4endl;

}
//-----------------------------------------------------------------------------
void Target::setY(G4double Y)
{

  Target_side_y=Y;
   G4cout<<"----> Target side Y is set to "<<G4BestUnit(Target_side_y,"Length")<<G4endl;

}
//-----------------------------------------------------------------------------
void Target::setZ(G4double Z)
{

  Target_thickness=Z;
  G4cout<<"----> Target thickness is set to "<<G4BestUnit(Target_thickness,"Length")<<G4endl;

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
  G4cout<<"----> Target side X is set to "<<G4BestUnit(2.*aTarget->GetXHalfLength(),"Length")<<G4endl;
  G4cout<<"----> Target side Y is set to "<<G4BestUnit(2.*aTarget->GetYHalfLength(),"Length")<<G4endl;  
  G4cout<<"----> Target thickness is set to "<<G4BestUnit(2.*aTarget->GetZHalfLength(),"Length")<<G4endl;
  G4cout<<"----> Number of simulation steps in the target is set to "<<NStep<<G4endl;
}
//---------------------------------------------------------------------
void Target::setMaterial(G4String materialName)
{

  TargetMaterialName = materialName;
  G4cout<<"----> Target material set to     "<<materialName<< G4endl;
         
}
//-------------------------------------------------------------------
void Target::setDensityScale(G4double scale)
{
  Target_density_scale = scale;
  G4cout<<"----> Scaling target density by     "<< scale << G4endl;
}
//-------------------------------------------------------------------
void Target::setTargetReactionDepth(G4double depth)
{
  //  G4cout<<"\n----> The depth is "<<G4BestUnit(depth,"Length")<< G4endl;;
  target_limits->SetUserMinRange(depth);
}
//-----------------------------------------------------------------------------
void Target::setPosition(G4double x, G4double y, G4double z)
{
  Pos->setX(x);
  Pos->setY(y);
  Pos->setZ(z);

  G4cout <<"----> Target position is set to "
	 << G4BestUnit(x,"Length") << ", "
    	 << G4BestUnit(y,"Length") << ", "
	 << G4BestUnit(z,"Length") << ", "
	 << G4endl;
}
//---------------------------------------------------------------------
void Target::setSourceFrame(G4String sF)
{
  if(Target_phys)
    delete Target_phys;
  
  sourceFrame = sF;

  if(sourceFrame == "eu152_Z2707"){

    frameMaterial = G4Material::GetMaterial("Al");
    //    frameThickness = 2.9*mm; // Used prior to 3/2012
    frameThickness = 1.2*mm;       // Dirk Weisshaar 3/4/2012
    frameInnerRadius = 3.8*cm/2.0; // Source data sheet
    frameOuterRadius = 5.4*cm/2.0; // Source data sheet
    tapeMaterial = G4Material::GetMaterial("G4_POLYETHYLENE"); // Source data sheet
    tapeThickness = 0.012*cm;      // Source data sheet

    euFrame = new G4Tubs("euFrame",frameInnerRadius,frameOuterRadius,frameThickness/2.,0.,360.*deg);
    euFrame_log = new G4LogicalVolume(euFrame,frameMaterial,"euFrame_log",0,0,0);
    euFrame_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),euFrame_log,"euFrame",expHall_log,false,0);

    euTape = new G4Tubs("euTape",0.,frameInnerRadius,tapeThickness/2.,0.,360.*deg);
    euTape_log = new G4LogicalVolume(euTape,tapeMaterial,"euTape_log",0,0,0);
    euTape_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),euTape_log,"euTape",expHall_log,false,0);

  } else if(sourceFrame == "cs137_E2879"){

    frameMaterial = G4Material::GetMaterial("Al");
    frameThickness = 0.7*mm;
    frameInnerRadius = 2.54*13./16.*cm/2.0;
    frameOuterRadius = 2.54*cm/2.0;
    tapeMaterial = G4Material::GetMaterial("G4_POLYETHYLENE");
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

  } else if(sourceFrame == "co56_2012"){

    frameMaterial = G4Material::GetMaterial("Fe");
    frameThickness = 0.05*mm;
    frameSide_x = 25.0*mm;
    frameSide_y = 25.0*mm;

    coFrame = new G4Box("coFrame", frameSide_x/2., frameSide_y/2., frameThickness/2.);
    coFrame_log = new G4LogicalVolume(coFrame,frameMaterial,"coFrame_log",0,0,0);
    coFrame_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos),coFrame_log,"coFrame",expHall_log,false,0);

  }

  G4cout<<"----> Including source frame: "<<sourceFrame<< G4endl;                 
}
//-------------------------------------------------------------------
void Target::BuildSled()
{

  sledMaterial = G4Material::GetMaterial("G10");

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
  *Pos0 += *Pos;
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
  *Pos0 += *Pos;
  sledRunner1_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), sledRunner_log, "sledRunner1_phys", expHall_log, false, 0);

  Pos0->set(0., -1.5*2.54*cm, 
	    sledBarLength/2.0 
	    + sledRunnerThickness/2.0 
	    + Target_thickness/2.0 
	    + sledFrameThickness/2.0
	    + 2.*tolerance);
  *Pos0 += *Pos;
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
//-------------------------------------------------------------------
#endif
