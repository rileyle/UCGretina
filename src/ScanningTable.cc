#ifdef SCANNING
#include "ScanningTable.hh"

ScanningTable::ScanningTable(G4LogicalVolume* experimentalHall_log,Materials* mat)
{
  CADModelPath = "../";

  materials=mat;
  expHall_log=experimentalHall_log;

  Pos0 = new G4ThreeVector(0., 0., 44.45*mm); // center flange @ (x,z) = (0,0)
  Pos2 = new G4ThreeVector(0., 0., -8.80*mm);

  includeCartFrame        = false;
  includeCloverCart       = false;
  includeSlitMount        = false;
  includeCollimator       = false;
  includeCollimatorInsert = false;
  includeShield           = false;
  includeCuTarget         = false;

  xShift  = 0.0*mm;
  yShift  = 0.0*mm;
  zShift  = 0.0*mm;
  cloverZ = 0.0*mm;
  //  cloverOffset = 343.88*mm; // Shield y position
  cloverOffset = 336.10*mm; // Mounts without shields

  collimatorRadius = 1.0*mm;

  material8020 = materials->FindMaterial("Al");
  materialCartBase = materials->FindMaterial("ssteel");
  materialCartTop = materials->FindMaterial("ssteel");
  materialCartBack = materials->FindMaterial("ssteel");
  materialSlitBrackets = materials->FindMaterial("Al");
  materialSlits = materials->FindMaterial("Hevimet");
  materialSlitAssembly = materials->FindMaterial("Al");
  materialTranslation = materials->FindMaterial("G10");
  materialTranslationAssembly = materials->FindMaterial("ssteel");
  materialCsCollimator = materials->FindMaterial("Hevimet");
  materialClover = materials->FindMaterial("HpGe");
  materialCloverShield = materials->FindMaterial("Al");
  materialRollers = materials->FindMaterial("G4_TEFLON");
  materialCuTarget = materials->FindMaterial("Cu");
}

ScanningTable::~ScanningTable()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* ScanningTable::Construct()
{

  //--- First the physical cart: 8020 frame, base, top and back panel ---------

  G4cout << "Loading STL files from directory:\n" << CADModelPath 
	 << "\nand building Tessellated solids ... " << G4endl;

  G4String CADFileName;

  if (includeCartFrame) {

    Pos3 = new G4ThreeVector(0., -360*mm, -424.45*mm);

    const int CartParts = 39;
    G4String CartPart[CartParts];
    G4Material* CartMaterial[CartParts];
    CartPart[0] = "Cart8020BackBrace1";
    CartMaterial[0] = material8020;
    CartPart[1] = "Cart8020BackBrace2";
    CartMaterial[1] = material8020;
    CartPart[2] = "Cart8020BackBrace3";
    CartMaterial[2] = material8020;
    CartPart[3] = "Cart8020BackBrace4";
    CartMaterial[3] = material8020;
    CartPart[4] = "Cart8020BackBrace5";
    CartMaterial[4] = material8020;
    CartPart[5] = "Cart8020BackBrace6";
    CartMaterial[5] = material8020;
    CartPart[6] = "Cart8020BackBrace7";
    CartMaterial[6] = material8020;
    CartPart[7] = "Cart8020BackBrace8";
    CartMaterial[7] = material8020;
    CartPart[8] = "Cart8020BottomBrace1";
    CartMaterial[8] = material8020;
    CartPart[9] = "Cart8020BottomBrace2";
    CartMaterial[9] = material8020;
    CartPart[10] = "Cart8020BottomBrace3";
    CartMaterial[10] = material8020;
    CartPart[11] = "Cart8020BottomBrace4";
    CartMaterial[11] = material8020;
    CartPart[12] = "Cart8020BottomBrace5";
    CartMaterial[12] = material8020;
    CartPart[13] = "Cart8020BottomBrace6";
    CartMaterial[13] = material8020;
    CartPart[14] = "Cart8020ExtBackLeft2";
    CartMaterial[14] = material8020;
    CartPart[15] = "Cart8020ExtBackMid2";
    CartMaterial[15] = material8020;
    CartPart[16] = "Cart8020ExtBackRight2";
    CartMaterial[16] = material8020;
    CartPart[17] = "Cart8020ExtBackTop2";
    CartMaterial[17] = material8020;
    CartPart[18] = "Cart8020ExtTopLeft2";
    CartMaterial[18] = material8020;
    CartPart[19] = "Cart8020ExtBottomLeft2";
    CartMaterial[19] = material8020;
    CartPart[20] = "Cart8020ExtBottomMid2";
    CartMaterial[20] = material8020;
    CartPart[21] = "Cart8020ExtBottomRight2";
    CartMaterial[21] = material8020;
    CartPart[22] = "Cart8020ExtFrontBottom2";
    CartMaterial[22] = material8020;
    CartPart[23] = "Cart8020ExtFrontLeft2";
    CartMaterial[23] = material8020;
    CartPart[24] = "Cart8020ExtFrontRight2";
    CartMaterial[24] = material8020;
    CartPart[25] = "Cart8020ExtTopRight2";
    CartMaterial[25] = material8020;
    CartPart[26] = "Cart8020FrontBrace1";
    CartMaterial[26] = material8020;
    CartPart[27] = "Cart8020FrontBrace2";
    CartMaterial[27] = material8020;
    CartPart[28] = "Cart8020LeftBrace1";
    CartMaterial[28] = material8020;
    CartPart[29] = "Cart8020LeftBrace2";
    CartMaterial[29] = material8020;
    CartPart[30] = "Cart8020LeftBrace3";
    CartMaterial[30] = material8020;
    CartPart[31] = "Cart8020LeftBrace4";
    CartMaterial[31] = material8020;
    CartPart[32] = "Cart8020RightBrace1";
    CartMaterial[32] = material8020;
    CartPart[33] = "Cart8020RightBrace2";
    CartMaterial[33] = material8020;
    CartPart[34] = "Cart8020RightBrace3";
    CartMaterial[34] = material8020;
    CartPart[35] = "Cart8020RightBrace4";
    CartMaterial[35] = material8020;
    CartPart[36] = "Cart8020ExtBottomBack2";
    CartMaterial[36] = material8020;
    CartPart[37] = "Cart8020BackBrace1";
    CartMaterial[37] = material8020;
    CartPart[38] = "Cart8020BackBrace2";
    CartMaterial[38] = material8020;

    G4Colour lpurple (0.5, 0.3, 1.0, 0.3);
    G4VisAttributes* VisCart = new G4VisAttributes(lpurple);
    VisCart->SetVisibility(true);
    VisCart->SetForceWireframe(true);

    for(int i=0;i<CartParts;i++){
      if(CartPart[i] != ""){
	CADFileName = CADModelPath + "/";
	CADFileName += CartPart[i];
	CADFileName += ".stl";
	CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
				    (char*)"STL");
	mesh->SetScale(mm);
	mesh->SetOffset(G4ThreeVector(0., 0., 0.));
     
	G4VSolid *Cart = mesh->TessellatedMesh();
	Cart_log = new G4LogicalVolume(Cart, CartMaterial[i], 
				       CartPart[i], 0, 0, 0);
	if(i<37)
	Cart_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
				      Cart_log,
				      CartPart[i],
				      expHall_log,false,0);
	else
	Cart_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos3),
				      Cart_log,
				      CartPart[i],
				      expHall_log,false,0);
	Cart_log->SetVisAttributes(VisCart);
      }
    }

    const int CartTopParts = 9;
    G4String CartTopPart[CartTopParts];
    G4Material* CartTopMaterial[CartTopParts];
    CartTopPart[0] = "CartTopBar";
    CartTopMaterial[0] = materialCartTop;
    CartTopPart[1] = "CartTopFlange";
    CartTopMaterial[1] = materialCartTop;
    CartTopPart[2] = "CartTopPlate";
    CartTopMaterial[2] = materialCartTop;
    CartTopPart[3] = "CartTopUBaseBottomLeft";
    CartTopMaterial[3] = materialCartTop;
    CartTopPart[4] = "CartTopUBaseBottomMid";
    CartTopMaterial[4] = materialCartTop;
    CartTopPart[5] = "CartTopUBaseBottomRight";
    CartTopMaterial[5] = materialCartTop;
    CartTopPart[6] = "CartTopUBaseUpperLeft";
    CartTopMaterial[6] = materialCartTop;
    CartTopPart[7] = "CartTopUBaseUpperMid";
    CartTopMaterial[7] = materialCartTop;
    CartTopPart[8] = "CartTopUBaseUpperRight";
    CartTopMaterial[8] = materialCartTop;

    for(int i=0; i < CartTopParts; i++){
      if(CartTopPart[i] != ""){
	CADFileName = CADModelPath + "/";
	CADFileName += CartTopPart[i];
	CADFileName += ".stl";
	CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
				    (char*)"STL");
	mesh->SetScale(mm);
	mesh->SetOffset(G4ThreeVector(0., 0., 0.));
     
	G4VSolid *CartTop = mesh->TessellatedMesh();
	
	CartTop_log = new G4LogicalVolume(CartTop,
					  CartTopMaterial[i], 
					  CartTopPart[i], 0, 0, 0);
	CartTop_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
					 CartTop_log,
					 CartTopPart[i],
					 expHall_log,false,0);
	CartTop_log->SetVisAttributes(VisCart);
      }
    }
  }
  
  //--- Now the slit assembly -------------------------------------------------

  // Slits on top
  // BotPos = new G4ThreeVector(0., 90.335*mm, -8.8*mm);
  // MidPos = new G4ThreeVector(0., 46.335*mm, -8.8*mm);
  // TopPos = new G4ThreeVector(0., -44.335*mm, -8.8*mm);

  // Slites on the bottom
  BotPos = new G4ThreeVector(0., 0., -8.8*mm);
  MidPos = new G4ThreeVector(0., 0., -8.8*mm);
  TopPos = new G4ThreeVector(0., 0., -8.8*mm);

  const int ZSlitParts = 10;
  G4String ZSlitPart[ZSlitParts];
  G4Material* ZSlitMaterial[ZSlitParts];
  ZSlitPart[0] = "ZSlitsBottomBracket";
  ZSlitMaterial[0] = materialSlitBrackets;
  ZSlitPart[1] = "ZSlitsLeftBracket";
  ZSlitMaterial[1] = materialSlitBrackets;
  ZSlitPart[2] = "ZSlitsRightBracket";
  ZSlitMaterial[2] = materialSlitBrackets;
  ZSlitPart[3] = "ZSlitsTopBracket";
  ZSlitMaterial[3] = materialSlitBrackets;
  ZSlitPart[4] = "ZSlitsLowerLeftWall";
  ZSlitMaterial[4] = materialSlits;
  ZSlitPart[5] = "ZSlitsLowerRightWall";
  ZSlitMaterial[5] = materialSlits;
  ZSlitPart[6] = "ZSlitsMidLeftWall";
  ZSlitMaterial[6] = materialSlits;
  ZSlitPart[7] = "ZSlitsMidRightWall";
  ZSlitMaterial[7] = materialSlits;
  ZSlitPart[8] = "ZSlitsUpperLeftWall";
  ZSlitMaterial[8] = materialSlits;
  ZSlitPart[9] = "ZSlitsUpperRightWall";
  ZSlitMaterial[9] = materialSlits;
  G4ThreeVector *UsePos;

  G4Colour lblue(0.0, 1.0, 1.0, 1.0);
  G4Colour lblue2(0.0, 1.0, 1.0, 0.3);

  G4VisAttributes* VisSlit = new G4VisAttributes(lblue);
  G4VisAttributes* VisSlit2 = new G4VisAttributes(lblue2);

  VisSlit->SetVisibility(true);
  VisSlit->SetForceSolid(false);
  VisSlit2->SetVisibility(true);
  VisSlit2->SetForceWireframe(true);
  
  for(int i=0;i<ZSlitParts;i++){
    if(ZSlitPart[i] != ""){
      CADFileName = CADModelPath + "/";
      CADFileName += ZSlitPart[i];
      CADFileName += ".stl";
      CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
				  (char*)"STL");
      mesh->SetScale(mm);
      mesh->SetOffset(G4ThreeVector(0., zShift, 0.));
     
      G4VSolid *Slits = mesh->TessellatedMesh();
      ZSlit_log = new G4LogicalVolume(Slits, ZSlitMaterial[i], 
				      ZSlitPart[i], 0, 0, 0);
	if(i == 4 || i == 5)
	  UsePos = BotPos;
	else if(i == 6 || i == 7)
	  UsePos = MidPos;
	else if(i == 8 || i == 9)
	  UsePos = TopPos;
	else
	  UsePos = Pos2;
      ZSlit_phys = new G4PVPlacement(G4Transform3D(NoRot, *UsePos),
				     ZSlit_log,
				     ZSlitPart[i],
				     expHall_log,false,0);
      if(i<4)
	ZSlit_log->SetVisAttributes(VisSlit2);
      else
	ZSlit_log->SetVisAttributes(VisSlit);
    }
  }

//-------------------------------------------------------------

  if(includeSlitMount){

    Pos1 = new G4ThreeVector(0., 0., 44.6*mm);

    const int SlitZAssemblyParts = 19;
    G4String SlitZAssemblyPart[SlitZAssemblyParts];
    G4Material* SlitZAssemblyMaterial[SlitZAssemblyParts];
    SlitZAssemblyPart[0] = "SlitZAssemblyBelt";
    SlitZAssemblyMaterial[0] = materialSlitAssembly;
    SlitZAssemblyPart[1] = "SlitZAssemblyCenterCylinder";
    SlitZAssemblyMaterial[1] = materialSlitAssembly;
    SlitZAssemblyPart[2] = "SlitZAssemblyCenterCylinderBase";
    SlitZAssemblyMaterial[2] = materialSlitAssembly;
    SlitZAssemblyPart[3] = "SlitZAssemblyCenterCylinderDisk";
    SlitZAssemblyMaterial[3] = materialSlitAssembly;
    SlitZAssemblyPart[4] = "SlitZAssemblyCenterCylinderPlug";
    SlitZAssemblyMaterial[4] = materialSlitAssembly;
    SlitZAssemblyPart[5] = "";
    SlitZAssemblyMaterial[5] = materialSlitAssembly;
    SlitZAssemblyPart[6] = "SlitZAssemblyLeftCrescentWithHoles";
    SlitZAssemblyMaterial[6] = materialSlitAssembly;
    SlitZAssemblyPart[7] = "SlitZAssemblyLeftBrace";
    SlitZAssemblyMaterial[7] = materialSlitAssembly;
    SlitZAssemblyPart[8] = "SlitZAssemblyLeftClamp";
    SlitZAssemblyMaterial[8] = materialSlitAssembly;
    SlitZAssemblyPart[9] = "SlitZAssemblyLeftPlate2";
    SlitZAssemblyMaterial[9] = materialSlitAssembly;
    SlitZAssemblyPart[10] = "SlitZAssemblyLeftPlate";
    SlitZAssemblyMaterial[10] = materialSlitAssembly;
    SlitZAssemblyPart[11] = "SlitZAssemblyPlateWithHoles";
    SlitZAssemblyMaterial[11] = materialSlitAssembly;
    SlitZAssemblyPart[12] = "";                          //  Get rid of this
    SlitZAssemblyMaterial[12] = materialSlitAssembly;    //  Get rid of this
    SlitZAssemblyPart[13] = "SlitZAssemblyRightBrace";
    SlitZAssemblyMaterial[13] = materialSlitAssembly;
    SlitZAssemblyPart[14] = "SlitZAssemblyRightClamp";
    SlitZAssemblyMaterial[14] = materialSlitAssembly;
    SlitZAssemblyPart[15] = "SlitZAssemblyRightCrescentWithHoles";
    SlitZAssemblyMaterial[15] = materialSlitAssembly;
    SlitZAssemblyPart[16] = "SlitZAssemblyRightPlate";
    SlitZAssemblyMaterial[16] = materialSlitAssembly;
    SlitZAssemblyPart[17] = "SlitZAssemblyTPlate";
    SlitZAssemblyMaterial[17] = materialSlitAssembly;
    SlitZAssemblyPart[18] = "SlitZAssemblyUPlate";
    SlitZAssemblyMaterial[18] = materialSlitAssembly;
    
    for(int i=0; i < SlitZAssemblyParts; i++){
      if(SlitZAssemblyPart[i] != ""){
	CADFileName = CADModelPath + "/";
	CADFileName += SlitZAssemblyPart[i];
	CADFileName += ".stl";
	CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
				    (char*)"STL");
	if((i != 0) && (   i < 5   || i == 6  || i == 8 || i == 9 
			|| i == 10 || i == 11 || i > 13)){
	mesh->SetScale(mm);
	mesh->SetOffset(G4ThreeVector(0., zShift, 0.));
	}
	else{
	mesh->SetScale(mm);
	mesh->SetOffset(G4ThreeVector(0., 0., 0.));
	}
     
	G4VSolid *SlitAssembly = mesh->TessellatedMesh();
	SlitZAssembly_log = new G4LogicalVolume(SlitAssembly, 
						SlitZAssemblyMaterial[i], 
						SlitZAssemblyPart[i], 0, 0, 0);
	if(i == 6 || i == 8 || i == 14 || i == 15){
	  SlitZAssembly_phys = new G4PVPlacement(G4Transform3D(NoRot, 
							       *Pos1),
						 SlitZAssembly_log,
						 SlitZAssemblyPart[i],
						 expHall_log,false,0);
	  SlitZAssembly_log->SetVisAttributes(VisSlit2);
	}
	else{
	  SlitZAssembly_phys = new G4PVPlacement(G4Transform3D(NoRot, 
							       *Pos0),
						 SlitZAssembly_log,
						 SlitZAssemblyPart[i],
						 expHall_log,false,0);
	  SlitZAssembly_log->SetVisAttributes(VisSlit2);
	}
  
	std::vector<G4TwoVector> polygon2;
	polygon2.push_back( G4TwoVector(0*mm,        0*mm) );
	polygon2.push_back( G4TwoVector(152.4*mm,    0*mm) );
	polygon2.push_back( G4TwoVector(152.4*mm,  -457.2*mm) );
	polygon2.push_back( G4TwoVector(133.35*mm, -457.2*mm) );
	polygon2.push_back( G4TwoVector(0*mm,      -19.05*mm) );

	std::vector<G4ExtrudedSolid::ZSection> zsections2;
	zsections2.push_back( G4ExtrudedSolid::ZSection(0.,       0., 1.) );
	zsections2.push_back( G4ExtrudedSolid::ZSection(19.05*mm, 0., 1.) );

	ZSlitAssemblyTriangle = new G4ExtrudedSolid("ZSlitAssemblyTriangle", 
						    polygon2,
						    zsections2);

	ZSlitAssemblyTriangle_log 
	  = new G4LogicalVolume(ZSlitAssemblyTriangle,
				SlitZAssemblyMaterial[i],
				"ZSlitAssemblyTriangle_log",
				0, 0, 0);
	ZSlitAssemblyTriangle_log->SetVisAttributes(VisSlit2);

	ZSlitAssemblyTriangle1Shift.setX(171.2*mm);
	ZSlitAssemblyTriangle1Shift.setY(397.5*mm + zShift);
	ZSlitAssemblyTriangle1Shift.setZ(-8.0*mm);
	ZSlitAssemblyTriangle1Pos = ZSlitAssemblyTriangle1Shift;

	ZSlitAssemblyTriangle2Shift.setX(-152.15*mm);
	ZSlitAssemblyTriangle2Shift.setY( 397.5*mm + zShift);
	ZSlitAssemblyTriangle2Shift.setZ( -8.0*mm);
	ZSlitAssemblyTriangle2Pos = ZSlitAssemblyTriangle2Shift;
  
	ZSlitAssemblyTriangleRot=G4RotationMatrix::IDENTITY;
	ZSlitAssemblyTriangleRot.rotateY(270.*deg);

	ZSlitAssemblyTriangle1_phys 
	  = new G4PVPlacement(G4Transform3D(ZSlitAssemblyTriangleRot,
					    ZSlitAssemblyTriangle1Pos),
			      ZSlitAssemblyTriangle_log, 
			      "ZSlitAssemblyTriangle1_phys",
			      expHall_log, false, 0);

	ZSlitAssemblyTriangle2_phys 
	  = new G4PVPlacement(G4Transform3D(ZSlitAssemblyTriangleRot,
					    ZSlitAssemblyTriangle2Pos),
			      ZSlitAssemblyTriangle_log, 
			      "ZSlitAssemblyTriangle2_phys",
			      expHall_log, false, 0);

      }
    }
  }

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

  if(includeCollimator){
    const int CollimatorParts = 2;
    G4String CollimatorPart[CollimatorParts];
    G4Material* CollimatorMaterial[CollimatorParts];
    CollimatorPart[0] = "CsCollimatorBase";
    CollimatorMaterial[0] = materialCsCollimator;
    CollimatorPart[1] = "CsCollimatorBody";
    CollimatorMaterial[1] = materialCsCollimator;
    //  CollimatorPart[2] = "CsCollimatorPlug";
    //  CollimatorMaterial[2] = materialCsCollimator;
    
    for(int i=0; i < CollimatorParts; i++){
      if(CollimatorPart[i] != ""){
	CADFileName = CADModelPath + "/";
	CADFileName += CollimatorPart[i];
	CADFileName += ".stl";
	CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
				    (char*)"STL");
	mesh->SetScale(mm);
	mesh->SetOffset(G4ThreeVector(xShift, 0., 37.05*mm + yShift));
     
	G4VSolid *Collimator = mesh->TessellatedMesh();

	Collimator_log = new G4LogicalVolume(Collimator,
					     CollimatorMaterial[i], 
					     CollimatorPart[i], 0, 0, 0);
	Collimator_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
					    Collimator_log,
					    CollimatorPart[i],
					    expHall_log,false,0);
	Collimator_log->SetVisAttributes(VisSlit2);
      }
    }

    if(includeCollimatorInsert){
      G4Tubs* CollimatorInsert = new G4Tubs("CollimatorInsert_vol",
					    collimatorRadius, 
					    7.5946*mm, 76.20*mm/2.0, 
					    0., 360.*deg);
      Collimator_log = new G4LogicalVolume(CollimatorInsert,
					   materialCsCollimator, 
					   "CollimatorInsert_log", 
					   0, 0, 0);
      // G4ThreeVector *collPos 
      // 	= new G4ThreeVector(Pos0->getX() + xShift, 
      // 			    Pos0->getY() + 249.301*mm, 
      // 			    Pos0->getZ() -  81.502*mm + yShift);
      G4ThreeVector *collIPos 
	= new G4ThreeVector(Pos0->getX() + xShift, 
			    Pos0->getY() + 249.301*mm, 
			    Pos0->getZ() - 44.45*mm + yShift);
      G4RotationMatrix collRot = G4RotationMatrix::IDENTITY;
      collRot.rotateX(90.*deg);
      Collimator_phys = new G4PVPlacement(G4Transform3D(collRot, *collIPos),
					  Collimator_log,
					  "CollimatorInsert",
					  expHall_log,false,0);
      Collimator_log->SetVisAttributes(VisSlit2);
    }
  }

  //--- Now the clover cart: base and elevator --------------------------------

  if (includeCloverCart) {

    // This is the mount for the BGO anti-compton shields (keep just un case they are ever used).
    
    // const int CloverAssemblyParts = 6;
    // G4String CloverAssemblyPart[CloverAssemblyParts];
    // G4Material* CloverAssemblyMaterial[CloverAssemblyParts];
    // CloverAssemblyPart[0] = "CloverAssemblyLeftBase";
    // CloverAssemblyMaterial[0] = materialCartBase;
    // CloverAssemblyPart[1] = "CloverAssemblyLeftLeft";
    // CloverAssemblyMaterial[1] = materialCartBase;
    // CloverAssemblyPart[2] = "CloverAssemblyLeftRight";
    // CloverAssemblyMaterial[2] = materialCartBase;
    // CloverAssemblyPart[3] = "CloverAssemblyRightBase";
    // CloverAssemblyMaterial[3] = materialCartBase;
    // CloverAssemblyPart[4] = "CloverAssemblyRightLeft";
    // CloverAssemblyMaterial[4] = materialCartBase;
    // CloverAssemblyPart[5] = "CloverAssemblyRightRight";
    // CloverAssemblyMaterial[5] = materialCartBase;

    // for(int i=0; i < CloverAssemblyParts; i++){
    //   if(CloverAssemblyPart[i] != ""){
    // 	CADFileName = CADModelPath + "/";
    // 	CADFileName += CloverAssemblyPart[i];
    // 	CADFileName += ".stl";
    // 	CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
    // 				    (char*)"STL");
    // 	mesh->SetScale(mm);
    // 	mesh->SetOffset(G4ThreeVector(xShift, 0., yShift));
	
    // 	G4VSolid *CloverAssembly = mesh->TessellatedMesh();
	
    // 	CloverAssembly_log = new G4LogicalVolume(CloverAssembly,
    // 						 CloverAssemblyMaterial[i], 
    // 						 CloverAssemblyPart[i], 0, 0, 0);
    // 	CloverAssembly_phys = new G4PVPlacement(G4Transform3D(NoRot, 
    // 							      *Pos0),
    // 						CloverAssembly_log,
    // 						CloverAssemblyPart[i],
    // 						expHall_log,false,0);
    // 	CloverAssembly_log->SetVisAttributes(VisSlit2);
    //   }
    // }

    const int CloverElevatorParts = 15;
    G4String CloverElevatorPart[CloverElevatorParts];
    G4Material* CloverElevatorMaterial[CloverElevatorParts];
    CloverElevatorPart[0] = "CloverElevatorBasePlate";
    CloverElevatorMaterial[0] = materialCartBase;
    CloverElevatorPart[1] = "CloverElevatorCubeLeftBack";
    CloverElevatorMaterial[1] = materialCartBase;
    CloverElevatorPart[2] = "CloverElevatorCubeLeftFront";
    CloverElevatorMaterial[2] = materialCartBase;
    CloverElevatorPart[3] = "CloverElevatorCubeRightBack";
    CloverElevatorMaterial[3] = materialCartBase;
    CloverElevatorPart[4] = "CloverElevatorCubeRightFront";
    CloverElevatorMaterial[4] = materialCartBase;
    CloverElevatorPart[5] = "CloverElevatorDogBone";
    CloverElevatorMaterial[5] = materialCartBase;
    CloverElevatorPart[6] = "CloverElevatorPlateLeft";
    CloverElevatorMaterial[6] = materialCartBase;
    CloverElevatorPart[7] = "CloverElevatorPlateRight";
    CloverElevatorMaterial[7] = materialCartBase;
    CloverElevatorPart[8] = "CloverElevatorRollerLeftBack";
    CloverElevatorMaterial[8] = materialRollers;
    CloverElevatorPart[9] = "CloverElevatorRollerLeftFront";
    CloverElevatorMaterial[9] = materialRollers;
    CloverElevatorPart[10] = "CloverElevatorRollerLeftMid";
    CloverElevatorMaterial[10] = materialRollers;
    CloverElevatorPart[11] = "CloverElevatorRollerRightBack";
    CloverElevatorMaterial[11] = materialRollers;
    CloverElevatorPart[12] = "CloverElevatorRollerRightFront";
    CloverElevatorMaterial[12] = materialRollers;
    CloverElevatorPart[13] = "CloverElevatorRollerRightMid";
    CloverElevatorMaterial[13] = materialRollers;
    CloverElevatorPart[14] = "CloverElevatorScrew";
    CloverElevatorMaterial[14] = materialCartBase;

    for(int i=0; i < CloverElevatorParts; i++){
      if(CloverElevatorPart[i] != ""){
	CADFileName = CADModelPath + "/";
	CADFileName += CloverElevatorPart[i];
	CADFileName += ".stl";
	CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
				    (char*)"STL");
	mesh->SetScale(mm);
	mesh->SetOffset(G4ThreeVector(0., cloverZ, 0.));
     
	G4VSolid *CloverElevator = mesh->TessellatedMesh();

	CloverElevator_log = new G4LogicalVolume(CloverElevator,
						 CloverElevatorMaterial[i], 
						 CloverElevatorPart[i], 0, 0, 0);
	CloverElevator_phys = new G4PVPlacement(G4Transform3D(NoRot, 
							      *Pos0),
						CloverElevator_log,
						CloverElevatorPart[i],
						expHall_log,false,0);
	CloverElevator_log->SetVisAttributes(VisSlit2);

	CloverMountBase = new G4Box("CloverMountBase", 
				    12.86*cm, 0.953*cm, 20.40*cm);

	CloverMountwall = new G4Box("CloverMountBase", 
				    12.86*cm, 8.017*cm, 0.953*cm);

	CloverMountCut = new G4Tubs("CloverMountCut",
				    0, 5.675*2*cm, 2.0*cm,
				    0*deg, 360.*deg);

	CloverMountSub 
	  = new G4SubtractionSolid("CloverMountSub", CloverMountwall, 
				   CloverMountCut, 
				   G4Transform3D(G4RotationMatrix(),
						 G4ThreeVector(0, 9.252*cm, 0)));

	std::vector<G4TwoVector> polygon;
	polygon.push_back( G4TwoVector(-12.86*cm, 0*cm) );
	polygon.push_back( G4TwoVector(-12.86*cm, 2.70*cm) );
	polygon.push_back( G4TwoVector(-8.10*cm, 15.527*cm) );
	polygon.push_back( G4TwoVector(-6.67*cm, 15.527*cm) );
	polygon.push_back( G4TwoVector(-6.67*cm,  9.02*cm) );
	polygon.push_back( G4TwoVector( 6.67*cm,  9.02*cm) );
	polygon.push_back( G4TwoVector( 6.66*cm, 15.527*cm) );
	polygon.push_back( G4TwoVector( 8.10*cm, 15.527*cm) );
	polygon.push_back( G4TwoVector(12.86*cm,  2.70*cm) );
	polygon.push_back( G4TwoVector(12.86*cm,  0*cm) );

	std::vector<G4ExtrudedSolid::ZSection> zsections;
	zsections.push_back( G4ExtrudedSolid::ZSection(-0.9525*cm, 0, 1) );
	zsections.push_back( G4ExtrudedSolid::ZSection( 0.9525*cm, 0, 1) );

	CloverMountExt = new G4ExtrudedSolid("CloverMountExt", polygon,
					     zsections);

	CloverMount 
	  = new G4UnionSolid("CloverMount", CloverMountBase, 
			     CloverMountSub,
			     G4Transform3D(G4RotationMatrix(),
					   G4ThreeVector(0, 8.017*cm, -19.45*cm)));

	CloverMount2 
	  = new G4UnionSolid("CloverMount2", CloverMount, CloverMountExt,
			     G4Transform3D(G4RotationMatrix(),
					   G4ThreeVector(0, 0*cm, 5.45*cm)));

	CloverMount_log = new G4LogicalVolume(CloverMount2,
					      materialCartBase,
					      "CloverMount_log",
					      0, 0, 0);
	CloverMount_log->SetVisAttributes(VisSlit2);

	CloverMountRot=G4RotationMatrix::IDENTITY;
	CloverMountRot.rotateY(20.*deg);

	CloverMountShift.setX(-19.7256*cm);
	CloverMountShift.setY(17.7*cm + cloverZ);
	CloverMountShift.setZ(-52.475*cm);
	CloverMountPos = CloverMountShift;

	CloverMount1Rot=G4RotationMatrix::IDENTITY;
	CloverMount1Rot.rotateY(-20.*deg);

	CloverMount1Shift.setX(19.7256*cm);
	CloverMount1Shift.setY(17.7*cm + cloverZ);
	CloverMount1Shift.setZ(-52.475*cm);
	CloverMount1Pos = CloverMount1Shift;

	CloverMount_phys = new G4PVPlacement(G4Transform3D(CloverMountRot, 
							   CloverMountPos),
					     CloverMount_log, "CloverMount_phys",
					     expHall_log, false, 0);

	CloverMount1_phys = new G4PVPlacement(G4Transform3D(CloverMount1Rot, 
							    CloverMount1Pos),
					      CloverMount_log, "CloverMount1_phys",
					      expHall_log, false, 0);

      }
    }
  }

 if (includeCuTarget) {
 
  CuTarget = new G4Box("CuTarget", 3.65/2*cm, 8.41375/2*cm, 3.49/2*cm);

  CuTarget_log = new G4LogicalVolume(CuTarget,
					      materialCuTarget,
					      "CuTarget_log",
					      0, 0, 0);

  CuTarget_log->SetVisAttributes(VisSlit2);

  CuTargetShift.setX(Pos0->getX() + xShift);
  CuTargetShift.setY(Pos0->getY() + 249.301*mm + 148.2625*mm);
  CuTargetShift.setZ(Pos0->getZ() -  44.45*mm + yShift);
  CuTargetPos = CuTargetShift;

  CuTarget_phys = new G4PVPlacement(G4Transform3D(NoRot, 
							    CuTargetPos),
					      CuTarget_log, "CuTarget_phys",
					      expHall_log, false, 0);
  }

    
  //--- Now the shields -- just for show right now ------------------------

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
    CloverRightShield_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos2),
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
    CloverLeftShield_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos2),
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
