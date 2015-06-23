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
  materialCsCollimator = materials->FindMaterial("Hevimet");
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

  G4String CADFileName = CADModelPath + "/CartTopUBaseBottomLeft.stl";
  CADMesh *meshCartTopUBBL = new CADMesh((char*)CADFileName.data(), 
					 (char*)"STL");
  meshCartTopUBBL->SetScale(mm);
  meshCartTopUBBL->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopUBBL = meshCartTopUBBL->TessellatedMesh();
  CartTopUBaseBottomLeft_log 
    = new G4LogicalVolume(CartTopUBBL, materialCartTop,
			  "CartTopUBaseBottomLeft_log", 0, 0, 0);
  CartTopUBaseBottomLeft_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
						  CartTopUBaseBottomLeft_log,
						  "CartTopUBaseBottomLeft", 
						  expHall_log, false, 0);
  CartTopUBaseBottomLeft_log->SetVisAttributes(VisCart);

  CADFileName = CADModelPath + "/CartTopUBaseBottomMid.stl";
  CADMesh *meshCartTopUBBM = new CADMesh((char*)CADFileName.data(), 
					 (char*)"STL");
  meshCartTopUBBM->SetScale(mm);
  meshCartTopUBBM->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopUBBM = meshCartTopUBBM->TessellatedMesh();
  CartTopUBaseBottomMid_log = new G4LogicalVolume(CartTopUBBM, 
						  materialCartTop,
						  "CartTopUBaseBottomMid_log", 
						  0, 0, 0);
  CartTopUBaseBottomMid_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
						 CartTopUBaseBottomMid_log,
						 "CartTopUBaseBottomMid", 
						 expHall_log, false, 0);
  CartTopUBaseBottomMid_log->SetVisAttributes(VisCart);

  CADFileName = CADModelPath + "/CartTopUBaseBottomRight.stl";
  CADMesh *meshCartTopUBBR = new CADMesh((char*)CADFileName.data(), 
					 (char*)"STL");
  meshCartTopUBBR->SetScale(mm);
  meshCartTopUBBR->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopUBBR = meshCartTopUBBR->TessellatedMesh();
  CartTopUBaseBottomRight_log 
    = new G4LogicalVolume(CartTopUBBR, materialCartTop,
			  "CartTopUBaseBottomRight_log", 
			  0, 0, 0);
  CartTopUBaseBottomRight_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
						   CartTopUBaseBottomRight_log,
						   "CartTopUBaseBottomRight", 
						   expHall_log, false, 0);
  CartTopUBaseBottomRight_log->SetVisAttributes(VisCart);

  CADFileName = CADModelPath + "/CartTopUBaseUpperLeft.stl";
  CADMesh *meshCartTopUBUL = new CADMesh((char*)CADFileName.data(), 
					 (char*)"STL");
  meshCartTopUBUL->SetScale(mm);
  meshCartTopUBUL->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopUBUL = meshCartTopUBUL->TessellatedMesh();
  CartTopUBaseUpperLeft_log = new G4LogicalVolume(CartTopUBUL, materialCartTop,
						  "CartTopUBaseUpperLeft_log", 
						  0, 0, 0);
  CartTopUBaseUpperLeft_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
						 CartTopUBaseUpperLeft_log,
						 "CartTopUBaseUpperLeft", 
						 expHall_log, false, 0);
  CartTopUBaseUpperLeft_log->SetVisAttributes(VisCart);

  CADFileName = CADModelPath + "/CartTopUBaseUpperMid.stl";
  CADMesh *meshCartTopUBUM = new CADMesh((char*)CADFileName.data(), 
					 (char*)"STL");
  meshCartTopUBUM->SetScale(mm);
  meshCartTopUBUM->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopUBUM = meshCartTopUBUM->TessellatedMesh();
  CartTopUBaseUpperMid_log = new G4LogicalVolume(CartTopUBUM, 
						 materialCartTop,
						 "CartTopUBaseUpperMid_log", 
						 0, 0, 0);
  CartTopUBaseUpperMid_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
						CartTopUBaseUpperMid_log,
						"CartTopUBaseUpperMid", 
						expHall_log, false, 0);
  CartTopUBaseUpperMid_log->SetVisAttributes(VisCart);

  CADFileName = CADModelPath + "/CartTopUBaseUpperRight.stl";
  CADMesh *meshCartTopUBUR = new CADMesh((char*)CADFileName.data(), 
					 (char*)"STL");
  meshCartTopUBUR->SetScale(mm);
  meshCartTopUBUR->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopUBUR = meshCartTopUBUR->TessellatedMesh();
  CartTopUBaseUpperRight_log = new G4LogicalVolume(CartTopUBUR, materialCartTop,
						   "CartTopUBaseUpperRight_log",
						   0, 0, 0);
  CartTopUBaseUpperRight_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
						  CartTopUBaseUpperRight_log,
						  "CartTopUBaseUpperRight", 
						  expHall_log, false, 0);
  CartTopUBaseUpperRight_log->SetVisAttributes(VisCart);

  CADFileName = CADModelPath + "/CartTopBar.stl";
  CADMesh *meshCartTopB = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshCartTopB->SetScale(mm);
  meshCartTopB->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopB = meshCartTopB->TessellatedMesh();
  CartTopBar_log = new G4LogicalVolume(CartTopB, materialCartTop,
				       "CartTopBar_log", 0,0,0);
  CartTopBar_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
				      CartTopBar_log,
				      "CartTopBar", expHall_log, false, 0);
  CartTopBar_log->SetVisAttributes(VisCart);

  CADFileName = CADModelPath + "/CartTopPlate.stl";
  CADMesh *meshCartTopP = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshCartTopP->SetScale(mm);
  meshCartTopP->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopP = meshCartTopP->TessellatedMesh();
  CartTopPlate_log = new G4LogicalVolume(CartTopP, materialCartTop,
					 "CartTopPlate_log", 0,0,0);
  CartTopPlate_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					CartTopPlate_log,
					"CartTopPlate", 
					expHall_log, false, 0);
  CartTopPlate_log->SetVisAttributes(VisCart);

  CADFileName = CADModelPath + "/CartTopFlange.stl";
  CADMesh *meshCartTopF = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshCartTopF->SetScale(mm);
  meshCartTopF->SetOffset(G4ThreeVector(0., 0., 0.));

  G4VSolid *CartTopF = meshCartTopF->TessellatedMesh();
  CartTopFlange_log = new G4LogicalVolume(CartTopF, materialCartTop,
					  "CartTopFlange_log", 0,0,0);
  CartTopFlange_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					 CartTopFlange_log,
					 "CartTopFlange", 
					 expHall_log, false, 0);
  CartTopFlange_log->SetVisAttributes(VisCart);

  // G4cout << "CartTop: center @" 
  // 	 << CartTop->GetExtent().GetExtentCenter() 
  // 	 << " xmin = " 
  // 	 << CartTop->GetExtent().GetXmin() 
  // 	 << " xmax = " 
  // 	 << CartTop->GetExtent().GetXmax() 
  // 	 << " ymin = " 
  // 	 << CartTop->GetExtent().GetYmin() 
  // 	 << " ymax = " 
  // 	 << CartTop->GetExtent().GetYmax() 
  // 	 << " zmin = " 
  // 	 << CartTop->GetExtent().GetZmin() 
  // 	 << " zmax = " 
  // 	 << CartTop->GetExtent().GetZmax() 
  // 	 << G4endl;

  if (includeCartFrame) {
    CADFileName = CADModelPath + "/Cart8020Frame.stl";
    CADMesh *mesh8020 = new CADMesh((char*)CADFileName.data(), 
				    (char*)"STL");
    mesh8020->SetScale(mm);
    mesh8020->SetOffset(G4ThreeVector(0., 0., 0.));

    CADFileName = CADModelPath + "/CartBaseOnly.stl";
    CADMesh *meshCartBase = new CADMesh((char*)CADFileName.data(), 
					(char*)"STL");
    meshCartBase->SetScale(mm);
    meshCartBase->SetOffset(G4ThreeVector(0., 0., 0.));

    CADFileName = CADModelPath + "/CartBack.stl";
    CADMesh *meshCartBack = new CADMesh((char*)CADFileName.data(), 
					(char*)"STL");
    meshCartBack->SetScale(mm);
    meshCartBack->SetOffset(G4ThreeVector(0., 0., 0.));
    
    CADFileName = CADModelPath + "/TurnTableBrake.stl";
    CADMesh *meshBrake = new CADMesh((char*)CADFileName.data(), 
				     (char*)"STL");
    meshBrake->SetScale(mm);
    meshBrake->SetOffset(G4ThreeVector(0., 0., 0.));
    
    G4VSolid *Cart8020 = mesh8020->TessellatedMesh();
    Cart8020_log = new G4LogicalVolume(Cart8020, material8020, 
				       "Cart8020_log", 0, 0, 0);
    Cart8020_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
				      Cart8020_log,
				      "Cart8020Frame", expHall_log, false, 0);
    Cart8020_log->SetVisAttributes(VisCart);
    
    G4VSolid *CartBase = meshCartBase->TessellatedMesh();
    CartBase_log = new G4LogicalVolume(CartBase, materialCartBase,
				       "CartBase_log", 0,0,0);
    CartBase_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
				      CartBase_log,
				      "CartBase", expHall_log, false, 0);
    CartBase_log->SetVisAttributes(VisCart);

    G4VSolid *CartBack = meshCartBack->TessellatedMesh();
    CartBack_log = new G4LogicalVolume(CartBack, materialCartBack,
				       "CartBack_log", 0, 0, 0);
    CartBack_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
				      CartBack_log,
				      "CartBack", expHall_log, false, 0);
    CartBack_log->SetVisAttributes(VisCart);
    
    G4VSolid *Brake = meshBrake->TessellatedMesh();
    Brake_log = new G4LogicalVolume(Brake, materialCartTop, "Brake_log", 
				    0, 0, 0);
    Brake_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
				   Brake_log, "BrakeFrame",
				   expHall_log, false, 0);
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

  CADFileName = CADModelPath + "/ZSlitsLeftBracket.stl";
  CADMesh *meshSlitsLB = new CADMesh((char*)CADFileName.data(), 
				     (char*)"STL");
  meshSlitsLB->SetScale(mm);
  meshSlitsLB->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsLB = meshSlitsLB->TessellatedMesh();
  ZSlitsLeftBracket_log = new G4LogicalVolume(SlitsLB, materialSlits,
					      "ZSlitsLeftBracket_log", 
					      0, 0, 0);
  ZSlitsLeftBracket_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					     ZSlitsLeftBracket_log,
					     "ZSlitsLeftBracket", 
					     expHall_log, false, 0);
  ZSlitsLeftBracket_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsUpperLeftWall.stl";
  CADMesh *meshSlitsULW = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshSlitsULW->SetScale(mm);
  meshSlitsULW->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsULW = meshSlitsULW->TessellatedMesh();
  ZSlitsUpperLeftWall_log = new G4LogicalVolume(SlitsULW, materialSlits,
						"ZSlitsUpperLeftWall", 
						0, 0, 0);
  ZSlitsUpperLeftWall_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					       ZSlitsUpperLeftWall_log,
					       "ZSlitsUpperLeftWall", 
					       expHall_log, false, 0);
  ZSlitsUpperLeftWall_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsMidLeftWall.stl";
  CADMesh *meshSlitsMLW = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshSlitsMLW->SetScale(mm);
  meshSlitsMLW->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsMLW = meshSlitsMLW->TessellatedMesh();
  ZSlitsMidLeftWall_log = new G4LogicalVolume(SlitsMLW, materialSlits,
					      "ZSlitsMidLeftWall", 
					      0, 0, 0);
  ZSlitsUpperLeftWall_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					       ZSlitsMidLeftWall_log,
					       "ZSlitsMidLeftWall", 
					       expHall_log, false, 0);
  ZSlitsMidLeftWall_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsLowerLeftWall.stl";
  CADMesh *meshSlitsLLW = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshSlitsLLW->SetScale(mm);
  meshSlitsLLW->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsLLW = meshSlitsLLW->TessellatedMesh();
  ZSlitsLowerLeftWall_log = new G4LogicalVolume(SlitsLLW, materialSlits,
						"ZSlitsLowerLeftWall", 
						0, 0, 0);
  ZSlitsLowerLeftWall_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					       ZSlitsLowerLeftWall_log,
					       "ZSlitsLowerLeftWall", 
					       expHall_log, false, 0);
  ZSlitsLowerLeftWall_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsUpperRightWall.stl";
  CADMesh *meshSlitsURW = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshSlitsURW->SetScale(mm);
  meshSlitsURW->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsURW = meshSlitsURW->TessellatedMesh();
  ZSlitsUpperRightWall_log = new G4LogicalVolume(SlitsURW, materialSlits,
						 "ZSlitsUpperRightWall", 
						 0, 0, 0);
  ZSlitsUpperRightWall_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
						ZSlitsUpperRightWall_log,
						"ZSlitsUpperRightWall", 
						expHall_log, false, 0);
  ZSlitsUpperRightWall_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsMidRightWall.stl";
  CADMesh *meshSlitsMRW = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshSlitsMRW->SetScale(mm);
  meshSlitsMRW->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsMRW = meshSlitsMRW->TessellatedMesh();
  ZSlitsMidRightWall_log = new G4LogicalVolume(SlitsMRW, materialSlits,
					       "ZSlitsMidRightWall", 
					       0, 0, 0);
  ZSlitsMidRightWall_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					      ZSlitsMidRightWall_log,
					      "ZSlitsMidRightWall", 
					      expHall_log, false, 0);
  ZSlitsMidRightWall_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsLowerRightWall.stl";
  CADMesh *meshSlitsLRW = new CADMesh((char*)CADFileName.data(), 
				      (char*)"STL");
  meshSlitsLRW->SetScale(mm);
  meshSlitsLRW->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsLRW = meshSlitsLRW->TessellatedMesh();
  ZSlitsLowerRightWall_log = new G4LogicalVolume(SlitsLRW, materialSlits,
						 "ZSlitsLowerRightWall", 
						 0, 0, 0);
  ZSlitsLowerRightWall_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
						ZSlitsLowerRightWall_log,
						"ZSlitsLowerRightWall", 
						expHall_log, false, 0);
  ZSlitsLowerRightWall_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsTopBracket.stl";
  CADMesh *meshSlitsTB = new CADMesh((char*)CADFileName.data(), 
				     (char*)"STL");
  meshSlitsTB->SetScale(mm);
  meshSlitsTB->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsTB = meshSlitsTB->TessellatedMesh();
  ZSlitsTopBracket_log = new G4LogicalVolume(SlitsTB, materialSlits,
					     "ZSlitsTopBracket", 
					     0, 0, 0);
  ZSlitsTopBracket_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					    ZSlitsTopBracket_log,
					    "ZSlitsTopBracket", 
					    expHall_log, false, 0);
  ZSlitsTopBracket_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsBottomBracket.stl";
  CADMesh *meshSlitsBB = new CADMesh((char*)CADFileName.data(), 
				     (char*)"STL");
  meshSlitsBB->SetScale(mm);
  meshSlitsBB->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsBB = meshSlitsBB->TessellatedMesh();
  ZSlitsBottomBracket_log = new G4LogicalVolume(SlitsBB, materialSlits,
						"ZSlitsBottomBracket", 
						0, 0, 0);
  ZSlitsBottomBracket_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					       ZSlitsBottomBracket_log,
					       "ZSlitsBottomBracket", 
					       expHall_log, false, 0);
  ZSlitsBottomBracket_log->SetVisAttributes(VisSlit);

  CADFileName = CADModelPath + "/ZSlitsRightBracket.stl";
  CADMesh *meshSlitsRB = new CADMesh((char*)CADFileName.data(), 
				     (char*)"STL");
  meshSlitsRB->SetScale(mm);
  meshSlitsRB->SetOffset(G4ThreeVector(0., zShift, 0.));
 
  G4VSolid *SlitsRB = meshSlitsRB->TessellatedMesh();
  ZSlitsRightBracket_log = new G4LogicalVolume(SlitsRB, materialSlits,
					       "ZSlitsRightBracket", 
					       0, 0, 0);
  ZSlitsRightBracket_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
					      ZSlitsRightBracket_log,
					      "ZSlitsRightBracket", 
					      expHall_log, false, 0);
  ZSlitsRightBracket_log->SetVisAttributes(VisSlit);

  // CADFileName = CADModelPath + "/SlitZAssemblyAll.stl";
  // CADMesh *meshSlitAssembly = new CADMesh((char*)CADFileName.data(), 
  // 					  (char*)"STL");
  // meshSlitAssembly->SetScale(mm);
  // meshSlitAssembly->SetOffset(G4ThreeVector(0., zShift, 0.));

  // G4VSolid *SlitAssembly = meshSlitAssembly->TessellatedMesh();
  // SlitAssembly_log = new G4LogicalVolume(SlitAssembly, materialSlitAssembly,
  // 					 "SlitAssembly_log", 0, 0, 0);
  // SlitAssembly_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
  // 					SlitAssembly_log,
  // 					"SlitAssembly", 
  // 					expHall_log, false, 0);
  // SlitAssembly_log->SetVisAttributes(VisSlit2);
  
  //--- And the Cs collimator assemblies --------------------------------------

  // G4Colour lgrey(0.7, 0.7, 0.7, 0.3);
  // G4VisAttributes* VisTranslation = new G4VisAttributes(lgrey);
  // VisTranslation->SetVisibility(true);
  // VisTranslation->SetForceWireframe(true);

  // CADFileName = CADModelPath + "/SourceTranslation1.stl";
  // CADMesh *meshTranslateX = new CADMesh((char*)CADFileName.data(), 
  // 					(char*)"STL");
  // meshTranslateX->SetScale(mm);
  // meshTranslateX->SetOffset(G4ThreeVector(xShift, 0., 0.));

  // CADFileName = CADModelPath + "/SourceTranslation2.stl";
  // CADMesh *meshTranslateY = new CADMesh((char*)CADFileName.data(), 
  // 					(char*)"STL");
  // meshTranslateY->SetScale(mm);
  // meshTranslateY->SetOffset(G4ThreeVector(0., 0., yShift));

  // CADFileName = CADModelPath + "/SourceTranslationAssembly.stl";
  // CADMesh *meshTranslate = new CADMesh((char*)CADFileName.data(), 
  // 				       (char*)"STL");
  // meshTranslate->SetScale(mm);
  // meshTranslate->SetOffset(G4ThreeVector(xShift, 0., yShift));

  // G4VSolid *TranslateX = meshTranslateX->TessellatedMesh();
  // TranslateX_log = new G4LogicalVolume(TranslateX, materialTranslation,
  // 				       "TranslateX_log", 0, 0, 0);
  // TranslateX_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
  // 				      TranslateX_log,
  // 				      "TranslationX", 
  // 				      expHall_log, false, 0);
  // TranslateX_log->SetVisAttributes(VisTranslation);

  // G4VSolid *TranslateY = meshTranslateY->TessellatedMesh();
  // TranslateY_log = new G4LogicalVolume(TranslateY, materialTranslation,
  // 				       "TranslateY_log", 0, 0, 0);
  // TranslateY_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
  // 				      TranslateY_log,
  // 				      "TranslationY", 
  // 				      expHall_log, false, 0);
  // TranslateY_log->SetVisAttributes(VisTranslation);
  
  // G4VSolid *Translate = meshTranslate->TessellatedMesh();
  // Translate_log = new G4LogicalVolume(Translate, materialTranslationAssembly,
  // 				      "Translate_log", 0, 0, 0);
  // Translate_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0), 
  // 				     Translate_log,
  // 				     "TranslationAssembly", 
  // 				     expHall_log, false, 0);
  // Translate_log->SetVisAttributes(VisTranslation);

  CADFileName = CADModelPath + "/CsCollimatorBase.stl";
  CADMesh *meshCsCollBase = new CADMesh((char*)CADFileName.data(), 
  					(char*)"STL");
  meshCsCollBase->SetScale(mm);
  meshCsCollBase->SetOffset(G4ThreeVector(xShift, 0., yShift));

  G4VSolid *CsCollimatorBase = meshCsCollBase->TessellatedMesh();
  CsCollimatorBase_log = new G4LogicalVolume(CsCollimatorBase, 
  					     materialCsCollimator,
  					     "CsCollimatorBase_log", 
					     0, 0, 0);
  CsCollimatorBase_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
  					    CsCollimatorBase_log, 
  					    "CsCollimatorBase",
  					    expHall_log, false, 0);
  CsCollimatorBase_log->SetVisAttributes(VisSlit2);

  CADFileName = CADModelPath + "/CsCollimatorBody.stl";
  CADMesh *meshCsCollBody = new CADMesh((char*)CADFileName.data(), 
  					(char*)"STL");
  meshCsCollBody->SetScale(mm);
  meshCsCollBody->SetOffset(G4ThreeVector(xShift, 0., yShift));

  G4VSolid *CsCollimatorBody = meshCsCollBody->TessellatedMesh();
  CsCollimatorBody_log = new G4LogicalVolume(CsCollimatorBody, 
  					     materialCsCollimator,
  					     "CsCollimatorBody_log", 
					     0, 0, 0);
  CsCollimatorBody_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
  					    CsCollimatorBody_log, 
  					    "CsCollimatorBody",
  					    expHall_log, false, 0);
  CsCollimatorBody_log->SetVisAttributes(VisSlit2);

  CADFileName = CADModelPath + "/CsCollimatorPlug.stl";
  CADMesh *meshCsCollPlug = new CADMesh((char*)CADFileName.data(), 
  					(char*)"STL");
  meshCsCollPlug->SetScale(mm);
  meshCsCollPlug->SetOffset(G4ThreeVector(xShift, 0., yShift));

  G4VSolid *CsCollimatorPlug = meshCsCollPlug->TessellatedMesh();
  CsCollimatorPlug_log = new G4LogicalVolume(CsCollimatorPlug, 
  					     materialCsCollimator,
  					     "CsCollimatorPlug_log", 
					     0, 0, 0);
  CsCollimatorPlug_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
  					    CsCollimatorPlug_log, 
  					    "CsCollimatorPlug",
  					    expHall_log, false, 0);
  CsCollimatorPlug_log->SetVisAttributes(VisSlit2);

  //--- Now the clover cart: base and elevator --------------------------------

  G4Colour lpurpleBlue (0.3, 0.2, 1.0, 0.3);
  G4VisAttributes* VisCloverCart = new G4VisAttributes(lpurpleBlue);
  VisCloverCart->SetVisibility(true);
  VisCloverCart->SetForceWireframe(true);

  if (includeCloverCart) {
    CADFileName = CADModelPath + "/CloverStand-Base.stl";
    CADMesh *meshCloverBase = new CADMesh((char*)CADFileName.data(), 
					  (char*)"STL");
    meshCloverBase->SetScale(mm);
    meshCloverBase->SetOffset(G4ThreeVector(0., 0., 0.));

    CADFileName = CADModelPath + "/CloverStand-Elevator.stl";
    CADMesh *meshCloverElevator = new CADMesh((char*)CADFileName.data(), 
					      (char*)"STL");
    meshCloverElevator->SetScale(mm);
    meshCloverElevator->SetOffset(G4ThreeVector(0., zShift, 0.));
    
    G4VSolid *CloverBase = meshCloverBase->TessellatedMesh();
    CloverBase_log = new G4LogicalVolume(CloverBase, materialCartBase,
					 "CloverBase_log", 0, 0, 0);
    CloverBase_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
					CloverBase_log, "CloverCartBase",
					expHall_log, false, 0);
    CloverBase_log->SetVisAttributes(VisCloverCart);
    
    G4VSolid *CloverElevator = meshCloverElevator->TessellatedMesh();
    CloverElevator_log = new G4LogicalVolume(CloverElevator, materialCartBase,
					     "CloverElevator_log", 
					     0, 0, 0);
    CloverElevator_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
					    CloverElevator_log, 
					    "CloverCartElevator",
					    expHall_log, false, 0);
    CloverElevator_log->SetVisAttributes(VisCloverCart);

    CADFileName = CADModelPath + "/CloverAssemblyRightBase.stl";
    CADMesh *meshCloverAssemblyRB = new CADMesh((char*)CADFileName.data(), 
						(char*)"STL");
    meshCloverAssemblyRB->SetScale(mm);
    meshCloverAssemblyRB->SetOffset(G4ThreeVector(0., zShift, 0.));
    
    G4VSolid *CloverAssemblyRB = meshCloverAssemblyRB->TessellatedMesh();
    CloverAssemblyRightBase_log 
      = new G4LogicalVolume(CloverAssemblyRB, materialCartBase, 
			    "CloverAssemblyRightBase_log", 0, 0, 0);
    CloverAssemblyRightBase_phys 
      = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
			  CloverAssemblyRightBase_log, 
			  "CloverAssemblyRightBase",
			  expHall_log, false, 0);
    CloverAssemblyRightBase_log->SetVisAttributes(VisCloverCart);

    CADFileName = CADModelPath + "/CloverAssemblyRightLeft.stl";
    CADMesh *meshCloverAssemblyRL = new CADMesh((char*)CADFileName.data(), 
						(char*)"STL");
    meshCloverAssemblyRL->SetScale(mm);
    meshCloverAssemblyRL->SetOffset(G4ThreeVector(0., zShift, 0.));
    
    G4VSolid *CloverAssemblyRL = meshCloverAssemblyRL->TessellatedMesh();
    CloverAssemblyRightLeft_log 
      = new G4LogicalVolume(CloverAssemblyRL, materialCartBase, 
			    "CloverAssemblyRightLeft_log", 0, 0, 0);
    CloverAssemblyRightLeft_phys 
      = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
			  CloverAssemblyRightLeft_log, 
			  "CloverAssemblyRightLeft",
			  expHall_log, false, 0);
    CloverAssemblyRightLeft_log->SetVisAttributes(VisCloverCart);

    CADFileName = CADModelPath + "/CloverAssemblyRightRight.stl";
    CADMesh *meshCloverAssemblyRR = new CADMesh((char*)CADFileName.data(), 
						(char*)"STL");
    meshCloverAssemblyRR->SetScale(mm);
    meshCloverAssemblyRR->SetOffset(G4ThreeVector(0., zShift, 0.));
    
    G4VSolid *CloverAssemblyRR = meshCloverAssemblyRR->TessellatedMesh();
    CloverAssemblyRightRight_log 
      = new G4LogicalVolume(CloverAssemblyRR, materialCartBase, 
			    "CloverAssemblyRightRight_log", 0, 0, 0);
    CloverAssemblyRightRight_phys 
      = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
			  CloverAssemblyRightRight_log, 
			  "CloverAssemblyRightRight",
			  expHall_log, false, 0);
    CloverAssemblyRightRight_log->SetVisAttributes(VisCloverCart);

    CADFileName = CADModelPath + "/CloverAssemblyLeftBase.stl";
    CADMesh *meshCloverAssemblyLB = new CADMesh((char*)CADFileName.data(), 
						(char*)"STL");
    meshCloverAssemblyLB->SetScale(mm);
    meshCloverAssemblyLB->SetOffset(G4ThreeVector(0., zShift, 0.));
    
    G4VSolid *CloverAssemblyLB = meshCloverAssemblyLB->TessellatedMesh();
    CloverAssemblyLeftBase_log 
      = new G4LogicalVolume(CloverAssemblyLB, materialCartBase, 
			    "CloverAssemblyLeftBase_log", 0, 0, 0);
    CloverAssemblyLeftBase_phys 
      = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
			  CloverAssemblyLeftBase_log, 
			  "CloverAssemblyLeftBase",
			  expHall_log, false, 0);
    CloverAssemblyLeftBase_log->SetVisAttributes(VisCloverCart);

    CADFileName = CADModelPath + "/CloverAssemblyLeftLeft.stl";
    CADMesh *meshCloverAssemblyLL = new CADMesh((char*)CADFileName.data(), 
						(char*)"STL");
    meshCloverAssemblyLL->SetScale(mm);
    meshCloverAssemblyLL->SetOffset(G4ThreeVector(0., zShift, 0.));
    
    G4VSolid *CloverAssemblyLL = meshCloverAssemblyLL->TessellatedMesh();
    CloverAssemblyLeftLeft_log 
      = new G4LogicalVolume(CloverAssemblyLL, materialCartBase, 
			    "CloverAssemblyLeftLeft_log", 0, 0, 0);
    CloverAssemblyLeftLeft_phys 
      = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
			  CloverAssemblyLeftLeft_log, 
			  "CloverAssemblyRightLeft",
			  expHall_log, false, 0);
    CloverAssemblyLeftLeft_log->SetVisAttributes(VisCloverCart);

    CADFileName = CADModelPath + "/CloverAssemblyLeftRight.stl";
    CADMesh *meshCloverAssemblyLR = new CADMesh((char*)CADFileName.data(), 
						(char*)"STL");
    meshCloverAssemblyLR->SetScale(mm);
    meshCloverAssemblyLR->SetOffset(G4ThreeVector(0., zShift, 0.));
    
    G4VSolid *CloverAssemblyLR = meshCloverAssemblyLR->TessellatedMesh();
    CloverAssemblyLeftRight_log 
      = new G4LogicalVolume(CloverAssemblyLR, materialCartBase, 
			    "CloverAssemblyLeftRight_log", 0, 0, 0);
    CloverAssemblyLeftRight_phys 
      = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
			  CloverAssemblyLeftRight_log, 
			  "CloverAssemblyLeftRight",
			  expHall_log, false, 0);
    CloverAssemblyLeftRight_log->SetVisAttributes(VisCloverCart);

  }
    
  //--- Now the clovers and shields -- just for show right now ----------------
  // CADFileName = CADModelPath + "/Clover-Right.stl";
  //  CADMesh *meshCloverRight = new CADMesh((char*)CADFileName.data(), 
  //					 (char*)"STL");
  // meshCloverRight->SetScale(mm);
  // meshCloverRight->SetOffset(G4ThreeVector(0., zShift, 0.));

  // G4VSolid *CloverRight = meshCloverRight->TessellatedMesh();
  // CloverRight_log = new G4LogicalVolume(CloverRight, materialClover, 
  // 					"CloverRight_log", 0, 0, 0);
  // CloverRight_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
  // 				       CloverRight_log, "RightClover",
  // 				       expHall_log, false, 0);

  // CADFileName = CADModelPath + "/Clover-Left.stl";
  // CADMesh *meshCloverLeft = new CADMesh((char*)CADFileName.data(), 
  // 					(char*)"STL");
  // meshCloverLeft->SetScale(mm);
  // meshCloverLeft->SetOffset(G4ThreeVector(0., zShift, 0.));

  // G4VSolid *CloverLeft = meshCloverLeft->TessellatedMesh();
  // CloverLeft_log = new G4LogicalVolume(CloverLeft, materialClover, 
  // 				       "CloverLeft_log", 0, 0, 0);
  // CloverLeft_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
  // 				       CloverLeft_log, "LeftClover",
  // 				       expHall_log, false, 0);

  if(includeShield){

    G4Colour lgreen (0.0, 1.0, 0.5, 0.3);
    G4VisAttributes* VisShield = new G4VisAttributes(lgreen);
    VisShield->SetVisibility(true);
    VisShield->SetForceWireframe(true);

    CADFileName = CADModelPath + "/Shield-Right.stl";
    CADMesh *meshCloverRightShield = new CADMesh((char*)CADFileName.data(), 
						 (char*)"STL");
    meshCloverRightShield->SetScale(mm);
    meshCloverRightShield->SetOffset(G4ThreeVector(0., zShift, 0.));

    G4VSolid *CloverRightShield = meshCloverRightShield->TessellatedMesh();
    CloverRightShield_log = new G4LogicalVolume(CloverRightShield, 
						materialCloverShield, 
						"CloverRightShield_log", 
						0, 0, 0);
    CloverRightShield_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
    					       CloverRightShield_log, 
    					       "RightCloverShield",
    					       expHall_log, false, 0);
    CloverRightShield_log->SetVisAttributes(VisShield);

    // G4cout << "CloverRightShield: center @" 
    // 	   << CloverRightShield->GetExtent().GetExtentCenter() 
    // 	   << " xmin = " 
    // 	   << CloverRightShield->GetExtent().GetXmin() 
    // 	   << " xmax = " 
    // 	   << CloverRightShield->GetExtent().GetXmax() 
    // 	   << " ymin = " 
    // 	   << CloverRightShield->GetExtent().GetYmin() 
    // 	   << " ymax = " 
    // 	   << CloverRightShield->GetExtent().GetYmax() 
    // 	   << " zmin = " 
    // 	   << CloverRightShield->GetExtent().GetZmin() 
    // 	   << " zmax = " 
    // 	   << CloverRightShield->GetExtent().GetZmax() 
    // 	   << G4endl;

    CADFileName = CADModelPath + "/Shield-Left.stl";
    CADMesh *meshCloverLeftShield = new CADMesh((char*)CADFileName.data(), 
						(char*)"STL");
    meshCloverLeftShield->SetScale(mm);
    meshCloverLeftShield->SetOffset(G4ThreeVector(0., zShift, 0.));

    G4VSolid *CloverLeftShield = meshCloverLeftShield->TessellatedMesh();
    CloverLeftShield_log = new G4LogicalVolume(CloverLeftShield, 
					       materialCloverShield, 
					       "CloverLeftShield_log", 
					       0, 0, 0);
    CloverLeftShield_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
    					      CloverLeftShield_log, 
    					      "LeftCloverShield",
    					      expHall_log, false, 0);
    CloverLeftShield_log->SetVisAttributes(VisShield);

  }

  G4cout << "Done." << G4endl;
  
  return Cart8020_phys;
}
//-----------------------------------------------------------------------------
#endif
