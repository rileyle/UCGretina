#ifdef SCANNING
#include "ScanningTable.hh"

ScanningTable::ScanningTable(G4LogicalVolume* experimentalHall_log,
			     Materials* mat)
{
  CADModelPath = "./cadModels/";

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

  slitWidth = 2.0*mm;

  material8020 = materials->FindMaterial("Al");
  materialCartBase = materials->FindMaterial("ssteel");
  materialCartTop = materials->FindMaterial("ssteel");
  materialCartBack = materials->FindMaterial("ssteel");
  materialSlitBrackets = materials->FindMaterial("Al");
  materialSlits = materials->FindMaterial("Hevimet");
  materialSlitAssembly = materials->FindMaterial("Al");
  materialSlitAssemblyCollimator = materials->FindMaterial("Cu");
  materialSlitAssemblyPlug = materials->FindMaterial("Hevimet");
  materialTranslation = materials->FindMaterial("G10");
  materialTranslationAssembly = materials->FindMaterial("ssteel");
  materialCsCollimator = materials->FindMaterial("Hevimet");
  materialCsCollimatorMount = materials->FindMaterial("Al");
  materialClover = materials->FindMaterial("Germanium");
  materialCloverShield = materials->FindMaterial("Al");
  materialRollers = materials->FindMaterial("G4_TEFLON");
  materialCuTarget = materials->FindMaterial("Cu");
}

ScanningTable::~ScanningTable()
{;}
//-----------------------------------------------------------------------------
void ScanningTable::Construct()
{

  //--- First the physical cart: 8020 frame, base, top and back panel ---------

  G4cout << "Constructing the LBL scanning table." << G4endl;

  G4cout << "   Loading STL files from directory:\n   " << CADModelPath 
	 << "\n   and building Tessellated solids ... " << G4endl;

  G4String CADFileName;

  if (includeCartFrame) {

    G4cout << "   ... cart frame ... " << G4endl;

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
	if(i<37) {
	  Cart_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
					Cart_log,
					CartPart[i],
					expHall_log,false,0);
	} else { // Front braces limiting the Clover stand vertically
	  Pos3 = new G4ThreeVector(0., -360*mm, -424.45*mm);
	  Cart_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos3),
					Cart_log,
					CartPart[i],
					expHall_log,false,0);
	}
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

  G4cout << "   ... slit assembly (slit width "
	 << slitWidth << " mm) ... " << G4endl;

  // Slits on the bottom
  BotPos = new G4ThreeVector(0., 0., -8.8*mm);
  MidPos = new G4ThreeVector(0., 0., -8.8*mm);
  TopPos = new G4ThreeVector(0., 0., -8.8*mm);

  // (The STL files have 2 mm slits.)
  G4ThreeVector *slitShift = new G4ThreeVector(0., slitWidth - 2.0*mm, 0.);
  (*MidPos) += (*slitShift);
  (*TopPos) += 2.0*(*slitShift);

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

      // if(i>3){
      // 	// This doubles the slit wall thickness with additional copies
      // 	G4ThreeVector* BP = new G4ThreeVector(0., 0., -8.8*mm + 16.22*mm);
      // 	G4ThreeVector* MP = new G4ThreeVector(0., 0., -8.8*mm + 16.22*mm);
      // 	G4ThreeVector* TP = new G4ThreeVector(0., 0., -8.8*mm + 16.22*mm);
      // 	if(i == 4 || i == 5)
      // 	  UsePos = BP;
      // 	else if(i == 6 || i == 7)
      // 	  UsePos = MP;
      // 	else if(i == 8 || i == 9)
      // 	  UsePos = TP;
      // 	// (The STL files have 2 mm slits.)
      // 	(*MP) += (*slitShift);
      // 	(*TP) += 2.0*(*slitShift);

      // 	G4VPhysicalVolume* ZS_phys 
      // 	  = new G4PVPlacement(G4Transform3D(NoRot, *UsePos),
      // 			      ZSlit_log,
      // 			      ZSlitPart[i],
      // 			      expHall_log,false,0);
      // 	// ... and this triples it.
      // 	BP->setZ(-8.8*mm + 2.0*16.22*mm);
      // 	MP->setZ(-8.8*mm + 2.0*16.22*mm);
      // 	TP->setZ(-8.8*mm + 2.0*16.22*mm);
      //        // (The STL files have 2 mm slits.)  
      //	(*MP) += (*slitShift);
      //	(*TP) += 2.0*(*slitShift);
      //
      // 	if(i == 4 || i == 5)
      // 	  UsePos = BP;
      // 	else if(i == 6 || i == 7)
      // 	  UsePos = MP;
      // 	else if(i == 8 || i == 9)
      // 	  UsePos = TP;
      // 	ZS_phys 
      // 	  = new G4PVPlacement(G4Transform3D(NoRot, *UsePos),
      // 			      ZSlit_log,
      // 			      ZSlitPart[i],
      // 			      expHall_log,false,0);
      //      }

      if(i<4)
	ZSlit_log->SetVisAttributes(VisSlit2);
      else
	ZSlit_log->SetVisAttributes(VisSlit);
    }
  }

//-------------------------------------------------------------

  if(includeSlitMount){

    G4cout << "   ... slit mount ... " << G4endl;

    Pos1 = new G4ThreeVector(0., 0., 44.6*mm);

    const int SlitZAssemblyParts = 18;
    G4String SlitZAssemblyPart[SlitZAssemblyParts];
    G4Material* SlitZAssemblyMaterial[SlitZAssemblyParts];
    SlitZAssemblyPart[0] = "SlitZAssemblyBelt";
    SlitZAssemblyMaterial[0] = materialSlitAssembly;
    SlitZAssemblyPart[1] = "SlitZAssemblyCenterCylinder";
    SlitZAssemblyMaterial[1] = materialSlitAssemblyCollimator;
    SlitZAssemblyPart[2] = "SlitZAssemblyCenterCylinderBase";
    SlitZAssemblyMaterial[2] = materialSlitAssemblyCollimator;
    SlitZAssemblyPart[3] = "SlitZAssemblyCenterCylinderDisk";
    SlitZAssemblyMaterial[3] = materialSlitAssemblyCollimator;
    SlitZAssemblyPart[4] = "SlitZAssemblyCenterCylinderPlug";
    SlitZAssemblyMaterial[4] = materialSlitAssemblyPlug;
    SlitZAssemblyPart[5] = "SlitZAssemblyLeftCrescentWithHoles";
    SlitZAssemblyMaterial[5] = materialSlitAssembly;
    SlitZAssemblyPart[6] = "SlitZAssemblyLeftBrace"; // Multiple solids
    SlitZAssemblyMaterial[6] = materialSlitAssembly;
    SlitZAssemblyPart[7] = "SlitZAssemblyLeftClamp";
    SlitZAssemblyMaterial[7] = materialSlitAssembly;
    SlitZAssemblyPart[8] = "SlitZAssemblyLeftPlate2";
    SlitZAssemblyMaterial[8] = materialSlitAssembly;
    //    SlitZAssemblyPart[9] = "SlitZAssemblyLeftPlate"; // overlaps TPlate
    SlitZAssemblyPart[9] = "";
    SlitZAssemblyMaterial[9] = materialSlitAssembly;
    SlitZAssemblyPart[10] = "SlitZAssemblyPlateWithHoles";
    SlitZAssemblyMaterial[10] = materialSlitAssembly;
    SlitZAssemblyPart[11] = "SlitZAssemblyRightBrace"; // Multiple solids
    SlitZAssemblyMaterial[11] = materialSlitAssembly;
    SlitZAssemblyPart[12] = "SlitZAssemblyRightClamp";
    SlitZAssemblyMaterial[12] = materialSlitAssembly;
    SlitZAssemblyPart[13] = "SlitZAssemblyRightCrescentWithHoles";
    SlitZAssemblyMaterial[13] = materialSlitAssembly;
    //    SlitZAssemblyPart[14] = "SlitZAssemblyRightPlate";
    SlitZAssemblyPart[14] = "";
    SlitZAssemblyMaterial[14] = materialSlitAssembly;
    SlitZAssemblyPart[15] = "SlitZAssemblyTPlate";
    SlitZAssemblyMaterial[15] = materialSlitAssembly;
    SlitZAssemblyPart[16] = "SlitZAssemblyUPlate";
    SlitZAssemblyMaterial[16] = materialSlitAssembly;
    SlitZAssemblyPart[17] = "CartBack";    
    SlitZAssemblyMaterial[17] = materialSlitAssembly;
    
    for(int i=0; i < SlitZAssemblyParts; i++){
      if(SlitZAssemblyPart[i] != ""){
	CADFileName = CADModelPath + "/";
	CADFileName += SlitZAssemblyPart[i];
	CADFileName += ".stl";
	CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
				    (char*)"STL");

        // Parts that translate
	if((i != 0) && (   i < 5  || i == 5  || i == 7 || i == 8 
			|| i == 9 || i == 10 || i > 11)){
	  mesh->SetScale(mm);
	  mesh->SetOffset(G4ThreeVector(0., zShift, 0.));
	} else { // and parts that don't
	  mesh->SetScale(mm);
	  mesh->SetOffset(G4ThreeVector(0., 0., 0.));
	}
     
	G4VSolid *SlitAssembly = mesh->TessellatedMesh();
	SlitZAssembly_log = new G4LogicalVolume(SlitAssembly, 
						SlitZAssemblyMaterial[i], 
						SlitZAssemblyPart[i], 
						0, 0, 0);
	if(i == 5 || i == 7 || i == 12 || i == 13){
	  SlitZAssembly_phys = new G4PVPlacement(G4Transform3D(NoRot, 
							       *Pos1),
						 SlitZAssembly_log,
						 SlitZAssemblyPart[i],
						 expHall_log,false,0);
	  SlitZAssembly_log->SetVisAttributes(VisSlit2);
	} else {
	  SlitZAssembly_phys = new G4PVPlacement(G4Transform3D(NoRot, 
							       *Pos0),
						 SlitZAssembly_log,
						 SlitZAssemblyPart[i],
						 expHall_log,false,0);
	  SlitZAssembly_log->SetVisAttributes(VisSlit2);
	}
      }
    }
    std::vector<G4TwoVector> polygon2;
    polygon2.push_back( G4TwoVector(0*mm,          0*mm) );
    polygon2.push_back( G4TwoVector(151.576*mm,    0*mm) );
    polygon2.push_back( G4TwoVector(151.576*mm, -457.2*mm) );
    polygon2.push_back( G4TwoVector(133.35*mm,  -457.2*mm) );
    polygon2.push_back( G4TwoVector(0*mm,        -19.05*mm) );

    std::vector<G4ExtrudedSolid::ZSection> zsections2;
    zsections2.push_back( G4ExtrudedSolid::ZSection(0.,       0., 1.) );
    zsections2.push_back( G4ExtrudedSolid::ZSection(19.05*mm, 0., 1.) );
	
    ZSlitAssemblyTriangle = new G4ExtrudedSolid("ZSlitAssemblyTriangle", 
						polygon2,
						zsections2);

    ZSlitAssemblyTriangle_log 
      = new G4LogicalVolume(ZSlitAssemblyTriangle,
			    materialSlitAssembly,
			    "ZSlitAssemblyTriangle_log",
			    0, 0, 0);
    ZSlitAssemblyTriangle_log->SetVisAttributes(VisSlit2);

    ZSlitAssemblyTriangle1Shift.setX(171.2*mm);
    ZSlitAssemblyTriangle1Shift.setY(397.38*mm + zShift);
    ZSlitAssemblyTriangle1Shift.setZ( -7.304*mm);
    ZSlitAssemblyTriangle1Pos = ZSlitAssemblyTriangle1Shift;

    ZSlitAssemblyTriangle2Shift.setX(-152.15*mm);
    ZSlitAssemblyTriangle2Shift.setY( 397.38*mm + zShift);
    ZSlitAssemblyTriangle2Shift.setZ(  -7.304*mm);
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

  //--- And the Cs collimator assemblies --------------------------------------

  if(includeCollimator){

    // Implementing the collimator as a stack of G4Tubs gives much better
    // navigation performance than tessellated solids and slightly
    // better than G4PolyCones.

    G4cout << "   ... collimator (x, y = " 
	   << xShift << ", " << yShift << " mm) ... " << G4endl;

    G4Tubs* CollimatorBase = new G4Tubs("CollimatorBase_vol",
					0., 2.0*25.4*mm,
					2.125*25.4*mm/2.0, 
					0., 360.*deg);

    CollimatorBase_log = new G4LogicalVolume(CollimatorBase,
					     materialCsCollimator, 
					     "CollimatorBase_log", 
					     0, 0, 0);

    G4ThreeVector *collBasePos 
      = new G4ThreeVector(Pos0->getX() + xShift, 
			  Pos0->getY() + 165.2095*mm, 
			  Pos0->getZ() - 44.45*mm + yShift);
    G4RotationMatrix collRot = G4RotationMatrix::IDENTITY;
    collRot.rotateX(90.*deg);
    CollimatorBase_phys = new G4PVPlacement(G4Transform3D(collRot, 
							  *collBasePos),
					    CollimatorBase_log,
					    "CollimatorBase",
					    expHall_log, false, 0);
    CollimatorBase_log->SetVisAttributes(VisSlit2);

    G4Tubs* CollimatorMid1 = new G4Tubs("CollimatorMid1_vol",
					0.75*25.4*mm, 2.0*25.4*mm,
					9.525*mm/2.0, 
					0., 360.*deg);

    CollimatorMid1_log = new G4LogicalVolume(CollimatorMid1,
					     materialCsCollimator, 
					     "CollimatorMid1_log", 
					     0, 0, 0);

    G4ThreeVector *collMid1Pos 
      = new G4ThreeVector(Pos0->getX() + xShift, 
			  Pos0->getY() + 196.9595*mm, 
			  Pos0->getZ() - 44.45*mm + yShift);
    CollimatorMid1_phys = new G4PVPlacement(G4Transform3D(collRot, 
    							  *collMid1Pos),
    					    CollimatorMid1_log,
    					    "CollimatorMid1",
    					    expHall_log, false, 0);
    CollimatorMid1_log->SetVisAttributes(VisSlit2);

    G4Tubs* CollimatorMid2 = new G4Tubs("CollimatorMid2_vol",
					1.005*25.4*mm, 2.0*25.4*mm,
					1.3*mm/2.0, 
					0., 360.*deg);

    CollimatorMid2_log = new G4LogicalVolume(CollimatorMid2,
					     materialCsCollimator, 
					     "CollimatorMid2_log", 
					     0, 0, 0);

    G4ThreeVector *collMid2Pos 
      = new G4ThreeVector(Pos0->getX() + xShift, 
			  Pos0->getY() + 202.372*mm, 
			  Pos0->getZ() - 44.45*mm + yShift);
    CollimatorMid2_phys = new G4PVPlacement(G4Transform3D(collRot, 
    							  *collMid2Pos),
    					    CollimatorMid2_log,
    					    "CollimatorMid2",
    					    expHall_log,false,0);
    CollimatorMid2_log->SetVisAttributes(VisSlit2);

    G4Tubs* CollimatorMid3 = new G4Tubs("CollimatorMid3_vol",
					0.125*25.4*mm, 2.0*25.4*mm,
					8.179*mm/2.0, 
					0., 360.*deg);

    CollimatorMid3_log = new G4LogicalVolume(CollimatorMid3,
					     materialCsCollimator, 
					     "CollimatorMid3_log", 
					     0, 0, 0);

    G4ThreeVector *collMid3Pos 
      = new G4ThreeVector(Pos0->getX() + xShift, 
			  Pos0->getY() + 207.1115*mm, 
			  Pos0->getZ() - 44.45*mm + yShift);
    CollimatorMid3_phys = new G4PVPlacement(G4Transform3D(collRot, 
    							  *collMid3Pos),
    					    CollimatorMid3_log,
    					    "CollimatorMid3",
    					    expHall_log,false,0);
    CollimatorMid3_log->SetVisAttributes(VisSlit2);

    G4Tubs* CollimatorBody = new G4Tubs("CollimatorBody_vol",
					0.31*25.4*mm/2.0, 2.0*25.4*mm,
					76.20*mm/2.0, 
					0., 360.*deg);

    CollimatorBody_log = new G4LogicalVolume(CollimatorBody,
					     materialCsCollimator, 
					     "CollimatorBody_log", 
					     0, 0, 0);

    G4ThreeVector *collBodyPos 
      = new G4ThreeVector(Pos0->getX() + xShift, 
			  Pos0->getY() + 249.301*mm, 
			  Pos0->getZ() - 44.45*mm + yShift);
    CollimatorBody_phys = new G4PVPlacement(G4Transform3D(collRot, 
    							  *collBodyPos),
    					    CollimatorBody_log,
    					    "CollimatorBoody",
    					    expHall_log, false, 0);
    CollimatorBody_log->SetVisAttributes(VisSlit2);

    if(includeCollimatorInsert){

      G4cout << "   ... collimator insert (" 
	     << collimatorRadius << " mm radius) ... " << G4endl;

      G4Tubs* CollimatorInsert = new G4Tubs("CollimatorInsert_vol",
					    collimatorRadius, 
					    0.3*25.4*mm/2.0, 76.20*mm/2.0, 
					    0., 360.*deg);
      Collimator_log = new G4LogicalVolume(CollimatorInsert,
					   materialCsCollimator, 
					   "CollimatorInsert_log", 
					   0, 0, 0);
      G4ThreeVector *collIPos 
	= new G4ThreeVector(Pos0->getX() + xShift, 
			    Pos0->getY() + 249.301*mm, 
			    Pos0->getZ() - 44.45*mm + yShift);
      collRot = G4RotationMatrix::IDENTITY;
      collRot.rotateX(90.*deg);
      Collimator_phys = new G4PVPlacement(G4Transform3D(collRot, *collIPos),
					  Collimator_log,
					  "CollimatorInsert",
					  expHall_log, false, 0);
      Collimator_log->SetVisAttributes(VisSlit2);
    }
  }

  if(includeCollimatorMount){

    G4cout << "   ... collimator mount ... " << G4endl;

    const int CollMountParts = 28;
    G4String CollMountPart[CollMountParts];
    G4Material* CollMountMaterial[CollMountParts];
    CollMountPart[0] = "SourceTranslationAssemblyBase";
    CollMountMaterial[0] = materialCsCollimatorMount;
    CollMountPart[1] = "SourceTranslationAssemblyClamp";
    CollMountMaterial[1] = materialCsCollimatorMount;
    CollMountPart[2] = "SourceTranslationAssemblyScrew";
    CollMountMaterial[2] = materialCsCollimatorMount;
    CollMountPart[3] = "SourceTranslationXStage";
    CollMountMaterial[3] = materialCsCollimatorMount;
    CollMountPart[4] = "SourceTranslationXBearing";
    CollMountMaterial[4] = materialCsCollimatorMount;
    CollMountPart[5] = "SourceTranslationXFrontPlate";
    CollMountMaterial[5] = materialCsCollimatorMount;
    CollMountPart[6] = "SourceTranslationXScrew";
    CollMountMaterial[6] = materialCsCollimatorMount;
    CollMountPart[7] = "SourceTranslationXBracket1";
    CollMountMaterial[7] = materialCsCollimatorMount;
    CollMountPart[8] = "SourceTranslationXMotorPlate";
    CollMountMaterial[8] = materialCsCollimatorMount;
    CollMountPart[9] = "SourceTranslationXBracket2";
    CollMountMaterial[9] = materialCsCollimatorMount;
    CollMountPart[10] = "SourceTranslationXMotor";
    CollMountMaterial[10] = materialCsCollimatorMount;
    CollMountPart[11] = "SourceTranslationXTrack";
    CollMountMaterial[11] = materialCsCollimatorMount;
    CollMountPart[12] = "SourceTranslationXEndPlate";
    CollMountMaterial[12] = materialCsCollimatorMount;
    CollMountPart[13] = "SourceTranslationXRing";
    CollMountMaterial[13] = materialCsCollimatorMount;
    CollMountPart[14] = "SourceTranslationYStage";
    CollMountMaterial[14] = materialCsCollimatorMount;
    CollMountPart[15] = "SourceTranslationYBasePlate";
    CollMountMaterial[15] = materialCsCollimatorMount;
    CollMountPart[16] = "SourceTranslationYBracket4";
    CollMountMaterial[16] = materialCsCollimatorMount;
    CollMountPart[17] = "SourceTranslationYRing";
    CollMountMaterial[17] = materialCsCollimatorMount;
    CollMountPart[18] = "SourceTranslationYBearing";
    CollMountMaterial[18] = materialCsCollimatorMount;
    CollMountPart[19] = "SourceTranslationYEndPlate";
    CollMountMaterial[19] = materialCsCollimatorMount;
    CollMountPart[20] = "SourceTranslationYScrew";
    CollMountMaterial[20] = materialCsCollimatorMount;
    CollMountPart[21] = "SourceTranslationYBracket1";
    CollMountMaterial[21] = materialCsCollimatorMount;
    CollMountPart[22] = "SourceTranslationYFrontPlate";
    CollMountMaterial[22] = materialCsCollimatorMount;
    CollMountPart[23] = "SourceTranslationYBracket2";
    CollMountMaterial[23] = materialCsCollimatorMount;
    CollMountPart[24] = "SourceTranslationYMotorPlate";
    CollMountMaterial[24] = materialCsCollimatorMount;
    CollMountPart[25] = "SourceTranslationYTrack";
    CollMountMaterial[25] = materialCsCollimatorMount;
    CollMountPart[26] = "SourceTranslationYBracket3";
    CollMountMaterial[26] = materialCsCollimatorMount;
    CollMountPart[27] = "SourceTranslationYMotor";
    CollMountMaterial[27] = materialCsCollimatorMount;

    for(int i=0; i<CollMountParts; i++){
      if(CollMountPart[i] != ""){
	CADFileName = CADModelPath + "/";
	CADFileName += CollMountPart[i];
	CADFileName += ".stl";
	CADMesh *mesh = new CADMesh((char*)CADFileName.data(),
				    (char*)"STL");
	mesh->SetScale(mm);
	if(i < 4)
	  mesh->SetOffset(G4ThreeVector(xShift, 0., 37.05*mm + yShift));
	else if(i < 15)
	  mesh->SetOffset(G4ThreeVector(0., 0., 37.05*mm + yShift));
	else if (i > 14)
	  mesh->SetOffset(G4ThreeVector(0., 0., 0.));

	G4VSolid *CollMount = mesh->TessellatedMesh();
	CollMount_log = new G4LogicalVolume(CollMount, CollMountMaterial[i], 
					    CollMountPart[i], 0, 0, 0);

	CollMount_phys = new G4PVPlacement(G4Transform3D(NoRot, *Pos0),
					   CollMount_log,
					   CollMountPart[i],
					   expHall_log, false, 0);
	CollMount_log->SetVisAttributes(VisSlit2);
      }
    }
  }

  //--- Now the clover cart: base and elevator --------------------------------

  if (includeCloverCart) {

    G4cout << "   ... Clover cart (z = " 
	   << zShift
	   << " mm) ... " << G4endl;

    // This is the mount for the BGO anti-compton shields.
    //      (keep just un case they are ever used)
    
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
    // 						 CloverAssemblyPart[i], 
    //                                           0, 0, 0);
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
						 CloverElevatorPart[i], 
						 0, 0, 0);
	CloverElevator_phys = new G4PVPlacement(G4Transform3D(NoRot, 
							      *Pos0),
						CloverElevator_log,
						CloverElevatorPart[i],
						expHall_log,false,0);
	CloverElevator_log->SetVisAttributes(VisSlit2);

      }
    }

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
					     G4ThreeVector(0, 
							   9.252*cm, 
							   0)));

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
				       G4ThreeVector(0, 
						     8.017*cm, 
						     -17.228*cm)));

    CloverMount2 
      = new G4UnionSolid("CloverMount2", CloverMount, CloverMountExt,
			 G4Transform3D(G4RotationMatrix(),
				       G4ThreeVector(0, 
						     0, 
						     2.69375*cm)));

    CloverMount_log = new G4LogicalVolume(CloverMount2,
					  materialCartBase,
					  "CloverMount_log",
					  0, 0, 0);
    CloverMount_log->SetVisAttributes(VisSlit2);

    CloverMountRot=G4RotationMatrix::IDENTITY;
    CloverMountRot.rotateY(20.*deg);

    CloverMountShift.setX(-20.4232*cm);
    CloverMountShift.setY(17.7*cm + cloverZ);
    CloverMountShift.setZ(-52.2643*cm);
    CloverMountPos = CloverMountShift;

    CloverMount1Rot=G4RotationMatrix::IDENTITY;
    CloverMount1Rot.rotateY(-20.*deg);
    
    CloverMount1Shift.setX(20.4232*cm);
    CloverMount1Shift.setY(17.7*cm + cloverZ);
    CloverMount1Shift.setZ(-52.2643*cm);
    CloverMount1Pos = CloverMount1Shift;

    CloverMount_phys = new G4PVPlacement(G4Transform3D(CloverMountRot, 
						       CloverMountPos),
					 CloverMount_log, 
					 "CloverMount_phys",
					 expHall_log, false, 0);

    CloverMount1_phys = new G4PVPlacement(G4Transform3D(CloverMount1Rot, 
							CloverMount1Pos),
					  CloverMount_log, 
					  "CloverMount1_phys",
					  expHall_log, false, 0);

  }

  if (includeCuTarget) {
    
    G4cout << "   ... Cu target ... " << G4endl;

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

    G4cout << "   ... clover anti-Compton shields ... " << G4endl;

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
  
  return;
}
//-----------------------------------------------------------------------------
#endif
