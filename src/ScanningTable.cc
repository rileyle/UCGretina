#ifdef SCANNING
#include "ScanningTable.hh"

ScanningTable::ScanningTable(G4LogicalVolume* experimentalHall_log,Materials* mat)
{
  CADModelPath = "../";

  materials=mat;
  expHall_log=experimentalHall_log;
  Pos0 = new G4ThreeVector(0., 0., 44.45*mm); // center flange @ (x,z) = (0,0)

  includeCartFrame  = false;
  includeCloverCart = false;
  includeShield     = false;

  xShift = 0.0*mm;
  yShift = 0.0*mm;
  zShift = 0.0*mm;
  
  material8020 = materials->FindMaterial("Al");
  materialCartBase = materials->FindMaterial("ssteel");
  materialCartTop = materials->FindMaterial("ssteel");
  materialCartBack = materials->FindMaterial("ssteel");
  materialSlits = materials->FindMaterial("Fe");
  materialSlitAssembly = materials->FindMaterial("ssteel");
  materialTranslation = materials->FindMaterial("G10");
  materialTranslationAssembly = materials->FindMaterial("ssteel");
  materialCsCollimator = materials->FindMaterial("Fe");
  materialClover = materials->FindMaterial("HpGe");
  materialCloverShield = materials->FindMaterial("Al");
}

ScanningTable::~ScanningTable()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* ScanningTable::Construct()
{

  //--- First the physical cart: 8020 frame, base, top and back panel ---------

  G4Colour lpurple (0.5, 0.3, 1.0, 0.3);
  G4VisAttributes* VisCart = new G4VisAttributes(lpurple);
  VisCart->SetVisibility(true);
  VisCart->SetForceWireframe(true);

  G4cout << "Loading STL files from directory:\n" << CADModelPath 
	 << "\nand building Tessellated solids ... " << G4endl;

  G4String CADFileName = CADModelPath + "/CartTopOnly.stl";
  //  G4String CADFileName = CADModelPath + "/CartFlangeOnly.stl";
  CADMesh *meshCartTop = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  meshCartTop->SetScale(mm);
  meshCartTop->SetOffset(G4ThreeVector(0.,0.,0.));

  G4VSolid *CartTop = meshCartTop->TessellatedMesh();
  CartTop_log = new G4LogicalVolume(CartTop, materialCartTop,
				     "CartTop_log", 0,0,0);
  CartTop_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), CartTop_log,
				    "CartTop", expHall_log, false, 0);
  CartTop_log->SetVisAttributes(VisCart);

  G4cout << "CartTop: center @" 
	 << CartTop->GetExtent().GetExtentCenter() 
	 << " xmin = " 
	 << CartTop->GetExtent().GetXmin() 
	 << " xmax = " 
	 << CartTop->GetExtent().GetXmax() 
	 << " ymin = " 
	 << CartTop->GetExtent().GetYmin() 
	 << " ymax = " 
	 << CartTop->GetExtent().GetYmax() 
	 << " zmin = " 
	 << CartTop->GetExtent().GetZmin() 
	 << " zmax = " 
	 << CartTop->GetExtent().GetZmax() 
	 << G4endl;

  if (includeCartFrame) {
    CADFileName = CADModelPath + "/Cart8020Frame.stl";
    CADMesh *mesh8020 = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    mesh8020->SetScale(mm);
    mesh8020->SetOffset(G4ThreeVector(0.,0.,0.));

    CADFileName = CADModelPath + "/CartBaseOnly.stl";
    CADMesh *meshCartBase = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshCartBase->SetScale(mm);
    meshCartBase->SetOffset(G4ThreeVector(0.,0.,0.));

    CADFileName = CADModelPath + "/CartBack.stl";
    CADMesh *meshCartBack = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshCartBack->SetScale(mm);
    meshCartBack->SetOffset(G4ThreeVector(0.,0.,0.));
    
    CADFileName = CADModelPath + "/TurnTableBrake.stl";
    CADMesh *meshBrake = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshBrake->SetScale(mm);
    meshBrake->SetOffset(G4ThreeVector(0.,0.,0.));
    
    G4VSolid *Cart8020 = mesh8020->TessellatedMesh();
    Cart8020_log = new G4LogicalVolume(Cart8020, material8020, "Cart8020_log", 0,0,0);
    Cart8020_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), Cart8020_log,
      "Cart8020Frame", expHall_log, false, 0);
    Cart8020_log->SetVisAttributes(VisCart);
    
    G4VSolid *CartBase = meshCartBase->TessellatedMesh();
    CartBase_log = new G4LogicalVolume(CartBase, materialCartBase,
      "CartBase_log", 0,0,0);
    CartBase_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), CartBase_log,
      "CartBase", expHall_log, false, 0);
    CartBase_log->SetVisAttributes(VisCart);

    G4VSolid *CartBack = meshCartBack->TessellatedMesh();
    CartBack_log = new G4LogicalVolume(CartBack, materialCartBack,
      "CartBack_log", 0,0,0);
    CartBack_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), CartBack_log,
      "CartBack", expHall_log, false, 0);
    CartBack_log->SetVisAttributes(VisCart);
    
    G4VSolid *Brake = meshBrake->TessellatedMesh();
    Brake_log = new G4LogicalVolume(Brake, materialCartTop,"Brake_log",0,0,0);
    Brake_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
      Brake_log,"BrakeFrame",expHall_log,false,0);
    Brake_log->SetVisAttributes(VisCart);
  }
  
  //--- Now the slit assembly -------------------------------------------------

  G4Colour lblue(0.0, 1.0, 1.0, 1.0);
  G4Colour lblue2(0.0, 1.0, 1.0, 0.3);
  G4VisAttributes* VisSlit = new G4VisAttributes(lblue);
  G4VisAttributes* VisSlit2 = new G4VisAttributes(lblue2);
  VisSlit->SetVisibility(true);
  VisSlit->SetForceSolid(false);
  VisSlit2->SetVisibility(true);
  VisSlit2->SetForceWireframe(true);

  CADFileName = CADModelPath + "/ZSlitsOnly.stl";
  CADMesh *meshSlits = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  meshSlits->SetScale(mm);
  meshSlits->SetOffset(G4ThreeVector(0.,zShift,0.));

  CADFileName = CADModelPath + "/SlitZAssemblyAll.stl";
  CADMesh *meshSlitAssembly = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  meshSlitAssembly->SetScale(mm);
  meshSlitAssembly->SetOffset(G4ThreeVector(0.,zShift,0.));

  G4VSolid *Slits = meshSlits->TessellatedMesh();
  Slits_log = new G4LogicalVolume(Slits, materialSlits,
				  "Slits_log", 0,0,0);
  Slits_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), Slits_log,
				 "Slits", expHall_log, false, 0);
  Slits_log->SetVisAttributes(VisSlit);

  G4VSolid *SlitAssembly = meshSlitAssembly->TessellatedMesh();
  SlitAssembly_log = new G4LogicalVolume(SlitAssembly, materialSlitAssembly,
					 "SlitAssembly_log", 0,0,0);
  SlitAssembly_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), SlitAssembly_log,
					"SlitAssembly", expHall_log, false, 0);
  SlitAssembly_log->SetVisAttributes(VisSlit2);
  
  //--- And the Cs collimator assemblies --------------------------------------

  G4Colour lgrey(0.7, 0.7, 0.7, 0.3);
  G4VisAttributes* VisTranslation = new G4VisAttributes(lgrey);
  VisTranslation->SetVisibility(true);
  VisTranslation->SetForceWireframe(true);

  CADFileName = CADModelPath + "/SourceTranslation1.stl";
  CADMesh *meshTranslateX = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  meshTranslateX->SetScale(mm);
  meshTranslateX->SetOffset(G4ThreeVector(xShift,0.,0.));

  CADFileName = CADModelPath + "/SourceTranslation2.stl";
  CADMesh *meshTranslateY = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  meshTranslateY->SetScale(mm);
  meshTranslateY->SetOffset(G4ThreeVector(0.,0.,yShift));

  CADFileName = CADModelPath + "/SourceTranslationAssembly.stl";
  CADMesh *meshTranslate = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  meshTranslate->SetScale(mm);
  meshTranslate->SetOffset(G4ThreeVector(xShift,0.,yShift));

  CADFileName = CADModelPath + "/CsCollimatorOnly.stl";
  CADMesh *meshCsColl = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  meshCsColl->SetScale(mm);
  meshCsColl->SetOffset(G4ThreeVector(xShift,0.,yShift));

  G4VSolid *TranslateX = meshTranslateX->TessellatedMesh();
  TranslateX_log = new G4LogicalVolume(TranslateX, materialTranslation,
				       "TranslateX_log", 0,0,0);
  TranslateX_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), TranslateX_log,
				      "TranslationX", expHall_log, false, 0);
  TranslateX_log->SetVisAttributes(VisTranslation);

  G4VSolid *TranslateY = meshTranslateY->TessellatedMesh();
  TranslateY_log = new G4LogicalVolume(TranslateY, materialTranslation,
				       "TranslateY_log", 0,0,0);
  TranslateY_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), TranslateY_log,
				      "TranslationY", expHall_log, false, 0);
  TranslateY_log->SetVisAttributes(VisTranslation);
  
  G4VSolid *Translate = meshTranslate->TessellatedMesh();
  Translate_log = new G4LogicalVolume(Translate, materialTranslationAssembly,
				      "Translate_log", 0,0,0);
  Translate_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0), Translate_log,
				     "TranslationAssembly", expHall_log, false, 0);
  Translate_log->SetVisAttributes(VisTranslation);

  G4VSolid *CsCollimator = meshCsColl->TessellatedMesh();
  CsCollimator_log = new G4LogicalVolume(CsCollimator, materialCsCollimator,
					 "CsCollimator_log", 0,0,0);
  CsCollimator_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
					CsCollimator_log, "CsCollimator",
					expHall_log, false, 0);
  CsCollimator_log->SetVisAttributes(VisSlit);


  G4cout << "CsCollimator: center @" 
	 << CsCollimator->GetExtent().GetExtentCenter() 
	 << " ymin = " 
	 << CsCollimator->GetExtent().GetYmin() 
	 << " ymax = " 
	 << CsCollimator->GetExtent().GetYmax() 
	 << " zmin = " 
	 << CsCollimator->GetExtent().GetZmin() 
	 << " zmax = " 
	 << CsCollimator->GetExtent().GetZmax() 
	 << G4endl;


  //--- Now the clover cart: base and elevator --------------------------------

  G4Colour lpurpleBlue (0.3, 0.2, 1.0, 0.3);
  G4VisAttributes* VisCloverCart = new G4VisAttributes(lpurpleBlue);
  VisCloverCart->SetVisibility(true);
  VisCloverCart->SetForceWireframe(true);

  if (includeCloverCart) {
    CADFileName = CADModelPath + "/CloverStand-Base.stl";
    CADMesh *meshCloverBase = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshCloverBase->SetScale(mm);
    meshCloverBase->SetOffset(G4ThreeVector(0.,0.,0.));

    CADFileName = CADModelPath + "/CloverStand-Elevator.stl";
    CADMesh *meshCloverElevator = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshCloverElevator->SetScale(mm);
    meshCloverElevator->SetOffset(G4ThreeVector(0.,zShift,0.));
    
    G4VSolid *CloverBase = meshCloverBase->TessellatedMesh();
    CloverBase_log = new G4LogicalVolume(CloverBase, materialCartBase,
					 "CloverBase_log", 0,0,0);
    CloverBase_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
					CloverBase_log, "CloverCartBase",
					expHall_log, false, 0);
    CloverBase_log->SetVisAttributes(VisCloverCart);
    
    G4VSolid *CloverElevator = meshCloverElevator->TessellatedMesh();
    CloverElevator_log = new G4LogicalVolume(CloverElevator, materialCartBase,
					     "CloverElevator_log", 0,0,0);
    CloverElevator_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
					    CloverElevator_log, "CloverCartElevator",
					    expHall_log, false, 0);
    CloverElevator_log->SetVisAttributes(VisCloverCart);

    CADFileName = CADModelPath + "/CloverAssembly-Right.stl";
    CADMesh *meshCloverAssemblyRight = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshCloverAssemblyRight->SetScale(mm);
    meshCloverAssemblyRight->SetOffset(G4ThreeVector(0.,zShift,0.));
    
    G4VSolid *CloverAssemblyRight = meshCloverAssemblyRight->TessellatedMesh();
    CloverAssemblyRight_log = new G4LogicalVolume(CloverAssemblyRight, materialCartBase, "CloverAssemblyRight_log", 0,0,0);
    CloverAssemblyRight_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
						 CloverAssemblyRight_log, "RightCloverAssembly",
						 expHall_log, false, 0);
    CloverAssemblyRight_log->SetVisAttributes(VisCloverCart);

    CADFileName = CADModelPath + "/CloverAssembly-Left.stl";
    CADMesh *meshCloverAssemblyLeft = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshCloverAssemblyLeft->SetScale(mm);
    meshCloverAssemblyLeft->SetOffset(G4ThreeVector(0.,zShift,0.));
    
    G4VSolid *CloverAssemblyLeft = meshCloverAssemblyLeft->TessellatedMesh();
    CloverAssemblyLeft_log = new G4LogicalVolume(CloverAssemblyLeft, materialCartBase, "CloverAssemblyLeft_log", 0,0,0);
    CloverAssemblyLeft_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
						CloverAssemblyLeft_log, "LeftCloverAssembly",
						expHall_log, false, 0);
    CloverAssemblyLeft_log->SetVisAttributes(VisCloverCart);
  }
    
  //--- Now the clovers and shields -- just for show right now ----------------
  // CADFileName = CADModelPath + "/Clover-Right.stl";
  // CADMesh *meshCloverRight = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  // meshCloverRight->SetScale(mm);
  // meshCloverRight->SetOffset(G4ThreeVector(0.,zShift,0.));

  // G4VSolid *CloverRight = meshCloverRight->TessellatedMesh();
  // CloverRight_log = new G4LogicalVolume(CloverRight, materialClover, "CloverRight_log", 0,0,0);
  // CloverRight_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
  // 				       CloverRight_log, "RightClover",
  // 				       expHall_log, false, 0);

  // CADFileName = CADModelPath + "/Clover-Left.stl";
  // CADMesh *meshCloverLeft = new CADMesh((char*)CADFileName.data(), (char*)"STL");
  // meshCloverLeft->SetScale(mm);
  // meshCloverLeft->SetOffset(G4ThreeVector(0.,zShift,0.));

  // G4VSolid *CloverLeft = meshCloverLeft->TessellatedMesh();
  // CloverLeft_log = new G4LogicalVolume(CloverLeft, materialClover, "CloverLeft_log", 0,0,0);
  // CloverLeft_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
  // 				       CloverLeft_log, "LeftClover",
  // 				       expHall_log, false, 0);

  if(includeShield){

    G4Colour lgreen (0.0, 1.0, 0.5, 0.3);
    G4VisAttributes* VisShield = new G4VisAttributes(lgreen);
    VisShield->SetVisibility(true);
    VisShield->SetForceWireframe(true);

    CADFileName = CADModelPath + "/Shield-Right.stl";
    CADMesh *meshCloverRightShield = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshCloverRightShield->SetScale(mm);
    meshCloverRightShield->SetOffset(G4ThreeVector(0.,zShift,0.));

    G4VSolid *CloverRightShield = meshCloverRightShield->TessellatedMesh();
    CloverRightShield_log = new G4LogicalVolume(CloverRightShield, materialCloverShield, "CloverRightShield_log", 0,0,0);
    CloverRightShield_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
					       CloverRightShield_log, "RightCloverShield",
					       expHall_log, false, 0);
    CloverRightShield_log->SetVisAttributes(VisShield);

    CADFileName = CADModelPath + "/Shield-Left.stl";
    CADMesh *meshCloverLeftShield = new CADMesh((char*)CADFileName.data(), (char*)"STL");
    meshCloverLeftShield->SetScale(mm);
    meshCloverLeftShield->SetOffset(G4ThreeVector(0.,zShift,0.));

    G4VSolid *CloverLeftShield = meshCloverLeftShield->TessellatedMesh();
    CloverLeftShield_log = new G4LogicalVolume(CloverLeftShield, materialCloverShield, "CloverLeftShield_log", 0,0,0);
    CloverLeftShield_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
					      CloverLeftShield_log, "LeftCloverShield",
  				       expHall_log, false, 0);
    CloverLeftShield_log->SetVisAttributes(VisShield);

  }

  G4cout << "Done." << G4endl;
  
  return Cart8020_phys;
}
//-----------------------------------------------------------------------------
#endif
