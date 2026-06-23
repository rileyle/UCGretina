#include "Greta_Shell.hh"

Greta_Shell::Greta_Shell()
{

  // Dimensions taken from fabricationprint_25j2106a-1.pdf
  Rmin      = 1022.0/2.0*mm;
  Rmax      = 1276.0/2.0*mm;

  smallPortRadius  = 7.0/2.0*2.54*cm;   // DETAIL C, SHEET2
  modulePortRadius = 12.008/2.0*2.54*cm; // DETAIL A, SHEET2
  northOffset      = 0.;
  southOffset      = 0.;

  // Small port positions
  PosSP[0] = G4ThreeVector(-492.891*mm, -145.253*mm, -256.924*mm);
  PosSP[1] = G4ThreeVector(-484.135*mm,  172.202*mm,  256.924*mm);
  PosSP[2] = G4ThreeVector(-313.380*mm, -407.226*mm,  256.924*mm);
  PosSP[3] = G4ThreeVector(-290.455*mm,  423.882*mm, -256.924*mm);
  PosSP[4] = G4ThreeVector( -14.168*mm, -513.653*mm, -256.924*mm);
  PosSP[5] = G4ThreeVector(  14.168*mm,  513.653*mm,  256.924*mm);
  PosSP[6] = G4ThreeVector( 290.455*mm, -423.882*mm,  256.924*mm);
  PosSP[7] = G4ThreeVector( 313.380*mm,  407.226*mm, -256.924*mm);
  PosSP[8] = G4ThreeVector( 484.135*mm, -172.202*mm, -256.924*mm);
  PosSP[9] = G4ThreeVector( 492.891*mm,  145.253*mm,  256.924*mm);

  // North: -1,  Split: 0,  South: 1, Omit: 2
  SmallPortStatus[0] =  2;
  SmallPortStatus[1] =  1;
  SmallPortStatus[2] =  1;
  SmallPortStatus[3] =  2;
  SmallPortStatus[4] =  0;
  SmallPortStatus[5] =  0;
  SmallPortStatus[6] =  2;
  SmallPortStatus[7] = -1;
  SmallPortStatus[8] = -1;
  SmallPortStatus[9] =  2;
  //Module LTriple positions
  /*
  MPosLTriple[0][0]=0*mm;
  MPosLTriple[0][1]=0*mm;
  MPosLTRiple[0][2]=0*mm;
  MPosLTriple[1][0]=0*mm;
  MPosLTriple[1][1]=0*mm;
  MPosLTriple[1][2]=0*mm;
  MPosLTriple[2][0]=0*mm;
  MPosLTriple[2][1]=0*mm;
  MPosLTriple[2][2]=0*mm;
  MPosLTriple[3][0]=0*mm;
  MPosLTriple[3][1]=0*mm;
  MPosLTriple[3][2]=0*mm;
  MPosLTriple[4][0]=0*mm;
  MPosLTriple[4][1]=0*mm;
  MPosLTriple[4][2]=0*mm;
  MPosLTriple[5][0]=0*mm;
  MPosLTriple[5][1]=0*mm;
  MPosLTriple[5][2]=0*mm;
  MPosLTriple[6][0]=0*mm;
  MPosLTriple[6][1]=0*mm;
  MPosLTriple[6][2]=0*mm;
  MPosLTriple[7][0]=0*mm;
  MPosLTriple[7][1]=0*mm;
  MPosLTriple[7][2]=0*mm;
  //Module LDouble positions
  MPosLDouble[0][0]=0*mm;
  MPosLDouble[0][1]=0*mm;
  MPosLDouble[0][2]=0*mm;
  MPosLDouble[1][0]=0*mm;
  MPosLDouble[1][1]=0*mm;
  MPosLDouble[1][2]=0*mm;
  MPosLDouble[2][0]=0*mm;
  MPosLDouble[2][1]=0*mm;
  MPosLDouble[2][2]=0*mm;
  MPosLDouble[3][0]=0*mm;
  MPosLDouble[3][1]=0*mm;
  MPosLDouble[3][2]=0*mm;
  MPosLDouble[4][0]=0*mm;
  MPosLDouble[4][1]=0*mm;
  MPosLDouble[4][2]=0*mm;
  MPosLDouble[5][0]=0*mm;
  MPosLDouble[5][1]=0*mm;
  MPosLDouble[5][2]=0*mm;
  MPosLDouble[6][0]=0*mm;
  MPosLDouble[6][1]=0*mm;
  MPosLDouble[6][2]=0*mm;
  MPosLDouble[7][0]=0*mm;
  MPosLDouble[7][1]=0*mm;
  MPosLDouble[7][2]=0*mm;
  //Module RTriple positions
  MPosRTriple[0][0]=0*mm;
  MPosRTriple[0][1]=0*mm;
  MPosRTRiple[0][2]=0*mm;
  MPosRTriple[1][0]=0*mm;
  MPosRTriple[1][1]=0*mm;
  MPosRTriple[1][2]=0*mm;
  MPosRTriple[2][0]=0*mm;
  MPosRTriple[2][1]=0*mm;
  MPosRTriple[2][2]=0*mm;
  MPosRTriple[3][0]=0*mm;
  MPosRTriple[3][1]=0*mm;
  MPosRTriple[3][2]=0*mm;
  MPosRTriple[4][0]=0*mm;
  MPosRTriple[4][1]=0*mm;
  MPosRTriple[4][2]=0*mm;
  MPosRTriple[5][0]=0*mm;
  MPosRTriple[5][1]=0*mm;
  MPosRTriple[5][2]=0*mm;
  MPosRTriple[6][0]=0*mm;
  MPosRTriple[6][1]=0*mm;
  MPosRTriple[6][2]=0*mm;
  MPosRTriple[7][0]=0*mm;
  MPosRTriple[7][1]=0*mm;
  MPosRTriple[7][2]=0*mm;
  //Module RDouble positions
  MPosRDouble[0][0]=0*mm;
  MPosRDouble[0][1]=0*mm;
  MPosRDouble[0][2]=0*mm;
  MPosRDouble[1][0]=0*mm;
  MPosRDouble[1][1]=0*mm;
  MPosRDouble[1][2]=0*mm;
  MPosRDouble[2][0]=0*mm;
  MPosRDouble[2][1]=0*mm;
  MPosRDouble[2][2]=0*mm;
  MPosRDouble[3][0]=0*mm;
  MPosRDouble[3][1]=0*mm;
  MPosRDouble[3][2]=0*mm;
  MPosRDouble[4][0]=0*mm;
  MPosRDouble[4][1]=0*mm;
  MPosRDouble[4][2]=0*mm;
  MPosRDouble[5][0]=0*mm;
  MPosRDouble[5][1]=0*mm;
  MPosRDouble[5][2]=0*mm;
  MPosRDouble[6][0]=0*mm;
  MPosRDouble[6][1]=0*mm;
  MPosRDouble[6][2]=0*mm;
  MPosRDouble[7][0]=0*mm;
  MPosRDouble[7][1]=0*mm;
  MPosRDouble[7][2]=0*mm;
  */
  //Module HexHole positions
  MPosHexHole[0][0]=-230*mm;
  MPosHexHole[0][1]=0*mm;
  MPosHexHole[0][2]=0*mm;
  MPosHexHole[1][0]=-230*cos(60*degree)*mm;
  MPosHexHole[1][1]=230*sin(60*degree)*mm;
  MPosHexHole[1][2]=0*mm;
  MPosHexHole[2][0]=230*cos(60*degree)*mm;
  MPosHexHole[2][1]=230*sin(60*degree)*mm;;
  MPosHexHole[2][2]=0*mm;
  MPosHexHole[3][0]=230*mm;
  MPosHexHole[3][1]=0*mm;
  MPosHexHole[3][2]=0*mm;
  //Module TripletHole positions
  MPosTripletHole[0][0]=-230*mm;
  MPosTripletHole[0][1]=0*mm;
  MPosTripletHole[0][2]=0*mm;
  MPosTripletHole[1][0]=-230*cos(60*degree)*mm;
  MPosTripletHole[1][1]=230*sin(60*degree)*mm;
  MPosTripletHole[1][2]=0*mm;
  MPosTripletHole[2][0]=230*cos(60*degree)*mm;
  MPosTripletHole[2][1]=230*sin(60*degree)*mm;
  MPosTripletHole[2][2]=0*mm;
  MPosTripletHole[3][0]=(1300/2*mm)*tan(31.717*degree);
  MPosTripletHole[3][1]=0*mm;
  MPosTripletHole[3][2]=0*mm;
  // Module Port Euler angles (relative to Slot 0)
  // Psi                                   Slot   Hemisphere
  // Theta
  // Phi
  ModuleEuler[0][0] =     0.000000*deg;   //  0    North, Left
  ModuleEuler[0][1] =    31.717473*deg;
  ModuleEuler[0][2] =   -54.000000*deg;
  ModuleEuler[1][0] =     0.000000*deg;   //  1    North, Left
  ModuleEuler[1][1] =    31.717473*deg;
  ModuleEuler[1][2] =    18.000000*deg;
  ModuleEuler[2][0] =     0.000000*deg;   //  2    Split
  ModuleEuler[2][1] =    31.717473*deg;
  ModuleEuler[2][2] =    90.000000*deg;
  ModuleEuler[3][0] =     0.000000*deg;   //  3    South, Right
  ModuleEuler[3][1] =    31.717473*deg;
  ModuleEuler[3][2] =   162.000000*deg;
  ModuleEuler[4][0] =     0.000000*deg;   //  4    South, Right
  ModuleEuler[4][1] =    31.717473*deg;
  ModuleEuler[4][2] =  -126.000000*deg;
  ModuleEuler[5][0] =    90.000000*deg;   //  5    North, Left
  ModuleEuler[5][1] =    58.282526*deg;
  ModuleEuler[5][2] =   -18.000000*deg;
  ModuleEuler[6][0] =    90.000000*deg;   //  6    North, Left
  ModuleEuler[6][1] =    58.282526*deg;
  ModuleEuler[6][2] =    54.000000*deg;
  ModuleEuler[7][0] =    90.000000*deg;   //  7    South, Right
  ModuleEuler[7][1] =    58.282526*deg;
  ModuleEuler[7][2] =   126.000000*deg;
  ModuleEuler[8][0] =    90.000000*deg;   //  8    South, Right
  ModuleEuler[8][1] =    58.282526*deg;
  ModuleEuler[8][2] =  -162.000000*deg;
  ModuleEuler[9][0] =    90.000000*deg;   //  9    Split
  ModuleEuler[9][1] =    58.282526*deg;
  ModuleEuler[9][2] =   -90.000000*deg;
  ModuleEuler[10][0] =  -31.750000*deg;   // 10    North, Left
  ModuleEuler[10][1] =   90.000000*deg;
  ModuleEuler[10][2] =  -72.000000*deg;
  ModuleEuler[11][0] =   31.750000*deg;   // 11    North, Left
  ModuleEuler[11][1] =   90.000000*deg;
  ModuleEuler[11][2] =  -36.000000*deg;
  ModuleEuler[12][0] =  -31.750000*deg;   // 12    North, Left
  ModuleEuler[12][1] =   90.000000*deg;
  ModuleEuler[12][2] =    0.000000*deg;
  ModuleEuler[13][0] =   31.750000*deg;   // 13    North, Left
  ModuleEuler[13][1] =   90.000000*deg;
  ModuleEuler[13][2] =   36.000000*deg;
  ModuleEuler[14][0] =  -31.750000*deg;   // 14    North, Left
  ModuleEuler[14][1] =   90.000000*deg;
  ModuleEuler[14][2] =   72.000000*deg;
  ModuleEuler[15][0] =   31.750000*deg;   // 15    South, Right
  ModuleEuler[15][1] =   90.000000*deg;
  ModuleEuler[15][2] =  108.000000*deg;
  ModuleEuler[16][0] =  -31.750000*deg;   // 16    South, Right
  ModuleEuler[16][1] =   90.000000*deg;
  ModuleEuler[16][2] =  144.000000*deg;
  ModuleEuler[17][0] =   31.750000*deg;   // 17    South, Right
  ModuleEuler[17][1] =   90.000000*deg;
  ModuleEuler[17][2] =  180.000000*deg;
  ModuleEuler[18][0] =  -31.750000*deg;   // 18    South, Right
  ModuleEuler[18][1] =   90.000000*deg;
  ModuleEuler[18][2] = -144.000000*deg;
  ModuleEuler[19][0] =   31.750000*deg;   // 19    South, Right
  ModuleEuler[19][1] =   90.000000*deg; 
  ModuleEuler[19][2] = -108.000000*deg;
  ModuleEuler[20][0] =  -90.000000*deg;   // 20    North, Left
  ModuleEuler[20][1] =  121.717475*deg;
  ModuleEuler[20][2] =  -54.000000*deg;
  ModuleEuler[21][0] =  -90.000000*deg;   // 21    North, Left
  ModuleEuler[21][1] =  121.717475*deg;
  ModuleEuler[21][2] =   18.000000*deg;
  ModuleEuler[22][0] =  -90.000000*deg;   // 22    Split, Left
  ModuleEuler[22][1] =  121.717475*deg;
  ModuleEuler[22][2] =   90.000000*deg;
  ModuleEuler[23][0] =  -90.000000*deg;   // 23    South, Right
  ModuleEuler[23][1] =  121.717475*deg;
  ModuleEuler[23][2] =  162.000000*deg;
  ModuleEuler[24][0] =  -90.000000*deg;   // 24    South, Right
  ModuleEuler[24][1] =  121.717475*deg;
  ModuleEuler[24][2] = -126.000000*deg;
  ModuleEuler[25][0] =  180.000000*deg;   // 25    North, Left
  ModuleEuler[25][1] =  148.282526*deg;
  ModuleEuler[25][2] =  -18.000000*deg;
  ModuleEuler[26][0] =  180.000000*deg;   // 26    North, Left
  ModuleEuler[26][1] =  148.282526*deg;
  ModuleEuler[26][2] =   54.000000*deg;
  ModuleEuler[27][0] =  180.000000*deg;   // 27    South, Right
  ModuleEuler[27][1] =  148.282526*deg;
  ModuleEuler[27][2] =  126.000000*deg;
  ModuleEuler[28][0] =  180.000000*deg;   // 28    South, Right
  ModuleEuler[28][1] =  148.282526*deg;
  ModuleEuler[28][2] = -162.000000*deg;
  ModuleEuler[29][0] =  180.000000*deg;   // 29    Split
  ModuleEuler[29][1] =  148.282526*deg;
  ModuleEuler[29][2] =  -90.000000*deg;

  // North: -1,  Split: 0,  South: 1
  ModulePortStatus[0]  = -1;
  ModulePortStatus[1]  = -1;
  ModulePortStatus[2]  = 0;
  ModulePortStatus[3]  = 1;
  ModulePortStatus[4]  = 1;
  ModulePortStatus[5]  = -1;
  ModulePortStatus[6]  = -1;
  ModulePortStatus[7]  = 1;
  ModulePortStatus[8]  = 1;
  ModulePortStatus[9]  = 0;
  ModulePortStatus[10]  = -1;
  ModulePortStatus[11]  = -1;
  ModulePortStatus[12]  = -1;
  ModulePortStatus[13]  = -1;
  ModulePortStatus[14]  = -1;
  ModulePortStatus[15]  = 1;
  ModulePortStatus[16]  = 1;
  ModulePortStatus[17]  = 1;
  ModulePortStatus[18]  = 1;
  ModulePortStatus[19]  = 1;
  ModulePortStatus[20]  = -1;
  ModulePortStatus[21]  = -1;
  ModulePortStatus[22]  = 0;
  ModulePortStatus[23]  = 1;
  ModulePortStatus[24]  = 1;
  ModulePortStatus[25]  = -1;
  ModulePortStatus[26]  = -1;
  ModulePortStatus[27]  = 1;
  ModulePortStatus[28]  = 1;
  ModulePortStatus[29]  = 0;

  matShell       = NULL;
  matShellName   = "Al";

  // This is the postion of the center of the front faces of a cluster before
  // it is placed according to its euler angles (the center of Slot 0).
  Pos0.setX(0.);
  Pos0.setY(0.);
  Pos0.setZ(180.0*mm);

  // Rescale to place the Tubs at a radial distance from the origin in the 
  // middle of the mounting shell. 
  Pos0.setMag((Rmin+Rmax)/2.*mm);

  Rot0 = G4RotationMatrix::IDENTITY;
  // This orients the cluster before it is placed according to its 
  // euler angles.
  Rot0.rotateY( Pos0.getTheta() );
  Rot0.rotateZ( Pos0.getPhi() );

}

Greta_Shell::~Greta_Shell()
{}

G4int Greta_Shell::FindMaterials()
{
  // search the material by its name
  G4Material* ptMaterial = G4Material::GetMaterial(matShellName);
  if (ptMaterial) {
    matShell = ptMaterial;
    G4cout << "\n ----> The GRETINA mounting shell material is "
          << matShell->GetName() << G4endl;
  }
  else {
    G4cout << " Could not find the material " << matShellName << G4endl;
    G4cout << " Could not build the GRETINA mounting shell! " << G4endl;
    return 1;
  } 
  return 0; 
}


void Greta_Shell::Placement(G4String status)
{
  if ( status != "Greta" && 
       status != "GretaLH" &&
       status != "Greta_North" && 
       status != "Greta_South" &&
       status != "GretaLH_North" && 
       status != "GretaLH_South"){
    G4cout << "Shell status " << status << " is not defined." << G4endl;
    return;
  }
  Test();
   
  /*
  if(status == "GretaLH" || status == "GretaLH_North")
    HalfShell("LH_north");

  if(status == "GretaLH" || status == "GretaLH_South")
    HalfShell("LH_south");
  */
   G4cout << "Constructed the " << status << " shell." << G4endl;
   G4cout << "  Shell radius: " << Rmin << " -- " << Rmax << G4endl;
   G4cout << "  Shell material: " << matShell->GetName() << G4endl;

}
//Forms and places a Greta Shell piece
void Greta_Shell::Test()
{
   G4RunManager* runManager = G4RunManager::GetRunManager();
   DetectorConstruction* theDetector = (DetectorConstruction*) runManager->GetUserDetectorConstruction();
   G4SubtractionSolid* shellL = Shell("LEFT");
   G4SubtractionSolid* shellR = Shell("RIGHT");
   //G4LogicalVolume* logicShell = new G4LogicalVolume(shell, matShell, "Shell_log", 0, 0, 0 );
   //new G4PVPlacement(0, G4ThreeVector(0,0,0), "MountingShell", logicShell, theDetector->HallPhys(), false, 0 );
   G4double halfheight = 1300/4*mm;
   G4double halfside = 230*mm;
   G4double innercut = 751*mm;
   G4ThreeVector NoShiftR = G4ThreeVector(0, 0, 0);
   G4RotationMatrix NoRotR = G4RotationMatrix::IDENTITY;;
   std::vector<G4TwoVector> polygon1(6);
   G4ThreeVector TarPosGlb = G4ThreeVector(0, 0, halfheight);
   G4RotationMatrix RotShellGlb = G4RotationMatrix::IDENTITY;
   polygon1[0] = G4TwoVector(MPosHexHole[0][0],MPosHexHole[0][1]);
   polygon1[1] = G4TwoVector(MPosHexHole[1][0],MPosHexHole[1][1]);
   polygon1[2] = G4TwoVector(MPosHexHole[2][0],MPosHexHole[2][1]);
   polygon1[3] = G4TwoVector(MPosHexHole[3][0],MPosHexHole[3][1]);
   polygon1[4] = G4TwoVector(halfside,-halfside);
   polygon1[5] = G4TwoVector(-halfside,-halfside);
   G4ExtrudedSolid* solidTarget1 = new G4ExtrudedSolid("solidTarget1",  polygon1, halfheight, G4TwoVector(0, 0), 0.00001, G4TwoVector(0, 0), 1);
   std::vector<G4TwoVector> polygon2(6);
   polygon2[0] = G4TwoVector(MPosTripletHole[0][0],MPosTripletHole[0][1]);
   polygon2[1] = G4TwoVector(MPosTripletHole[1][0],MPosTripletHole[1][1]);
   polygon2[2] = G4TwoVector(MPosTripletHole[2][0],MPosTripletHole[2][1]);
   polygon2[3] = G4TwoVector(MPosTripletHole[3][0],MPosTripletHole[3][1]);
   polygon2[4] = G4TwoVector((1300/2*mm)*tan(31.717*degree),-halfside);
   polygon2[5] = G4TwoVector(-halfside,-halfside);
   G4ExtrudedSolid* solidTarget2 = new G4ExtrudedSolid("solidTarget2",  polygon2, halfheight, G4TwoVector(0, 0), 0.00001, G4TwoVector(0, 0), 1);
 
   G4ThreeVector TarPos10 = G4ThreeVector(0, 0, halfheight);
   TarPos10.rotateZ(ModuleEuler[9][0]);
   TarPos10.rotateY(ModuleEuler[9][1]);
   TarPos10.rotateZ(ModuleEuler[9][2]);
   G4RotationMatrix RotShell10 = G4RotationMatrix::IDENTITY;
   RotShell10.rotateZ( 0*degree );
   RotShell10.rotateY( TarPos10.getTheta() );
   RotShell10.rotateZ( TarPos10.getPhi() );
   G4SubtractionSolid* cutout10 = new G4SubtractionSolid("cutout10", shellL, solidTarget1, G4Transform3D(RotShell10, TarPos10));

   // Goes with the Right "hemisphere"
   G4IntersectionSolid* bump10 = new G4IntersectionSolid("bump10",shellL, solidTarget1, G4Transform3D(RotShell10, TarPos10));
   G4LogicalVolume* logicbump10 = new G4LogicalVolume(bump10, matShell, "Shell_log", 0, 0, 0 );

   G4ThreeVector TarPos30 = G4ThreeVector(0, 0, halfheight);
   TarPos30.rotateZ(ModuleEuler[29][0]);
   TarPos30.rotateY(ModuleEuler[29][1]);
   TarPos30.rotateZ(ModuleEuler[29][2]);
   G4RotationMatrix RotShell30 = G4RotationMatrix::IDENTITY;
   RotShell30.rotateZ(0*degree );
   RotShell30.rotateY( TarPos30.getTheta() );
   RotShell30.rotateZ( TarPos30.getPhi() );
   // Left hemisphere with both cutouts but not the bumpouts
   G4SubtractionSolid* cutout30 = new G4SubtractionSolid("cutout30", cutout10, solidTarget2, G4Transform3D(RotShell30, TarPos30));

   // Goes with the Right "hemisphere"
   G4IntersectionSolid* bump30 = new G4IntersectionSolid("bump30",cutout10, solidTarget2, G4Transform3D(RotShell30, TarPos30));
   G4LogicalVolume* logicbump30 = new G4LogicalVolume(bump30, matShell, "Shell_log", 0, 0, 0 );
   G4LogicalVolume* logicCutout30 = new G4LogicalVolume(cutout30, matShell, "Shell_log", 0, 0, 0 );
   
   G4ThreeVector TarPos23 = G4ThreeVector(0, 0, halfheight);
   TarPos23.rotateZ(ModuleEuler[22][0]);
   TarPos23.rotateY(ModuleEuler[22][1]);
   TarPos23.rotateZ(ModuleEuler[22][2]);
   G4RotationMatrix RotShell23 = G4RotationMatrix::IDENTITY;
   RotShell23.rotateZ( 0*degree );
   RotShell23.rotateY( TarPos23.getTheta() );
   RotShell23.rotateZ( TarPos23.getPhi() );
   // Right hemisphere with the hole 23 cutout
   G4SubtractionSolid* cutout23 = new G4SubtractionSolid("cutout23", shellR, solidTarget1, G4Transform3D(RotShell23, TarPos23));

   // Goes with the Left hemisphere.
   G4IntersectionSolid* bump23 = new G4IntersectionSolid("bump23", shellR, solidTarget1, G4Transform3D(RotShell23, TarPos23));
   G4LogicalVolume* logicbump23 = new G4LogicalVolume(bump23, matShell, "Shell_log", 0, 0, 0 );
   
   G4ThreeVector TarPos3 = G4ThreeVector(0, 0, halfheight);
   TarPos3.rotateZ(ModuleEuler[2][0]);
   TarPos3.rotateY(ModuleEuler[2][1]);
   TarPos3.rotateZ(ModuleEuler[2][2]);
   G4RotationMatrix RotShell3 = G4RotationMatrix::IDENTITY;
   RotShell3.rotateZ( 0*degree );
   RotShell3.rotateY( TarPos3.getTheta() );
   RotShell3.rotateZ( TarPos3.getPhi() );
   // Right hemisphere with both cutouts but not the bumpouts
   G4SubtractionSolid* cutout3 = new G4SubtractionSolid("cutout3", cutout23, solidTarget2, G4Transform3D(RotShell3, TarPos3));
   G4LogicalVolume* logicCutout3 = new G4LogicalVolume(cutout3, matShell, "Shell_log", 0, 0, 0 );
   
   // Goes with the Left hemisphere
   G4IntersectionSolid* bump3 = new G4IntersectionSolid("bump3", cutout23, solidTarget2, G4Transform3D(RotShell3, TarPos3));
   G4LogicalVolume* logicbump3 = new G4LogicalVolume(bump3, matShell, "Shell_log", 0, 0, 0 );
   
   G4ThreeVector NoShiftL = G4ThreeVector(0, 0, 0);
   G4RotationMatrix NoRotL = G4RotationMatrix::IDENTITY;
   
   std::vector<G4TwoVector> polygon3(5);
   polygon3[0] = G4TwoVector(innercut*sin(ModuleEuler[2][1])*cos(ModuleEuler[2][2]),innercut*sin(ModuleEuler[2][1])*sin(ModuleEuler[2][2]));
   polygon3[1] = G4TwoVector(innercut*sin(ModuleEuler[1][1])*cos(ModuleEuler[1][2]),innercut*sin(ModuleEuler[1][1])*sin(ModuleEuler[1][2]));
   polygon3[2] = G4TwoVector(innercut*sin(ModuleEuler[0][1])*cos(ModuleEuler[0][2]),innercut*sin(ModuleEuler[0][1])*sin(ModuleEuler[0][2]));
   polygon3[3] = G4TwoVector(innercut*sin(ModuleEuler[4][1])*cos(ModuleEuler[4][2]),innercut*sin(ModuleEuler[4][1])*sin(ModuleEuler[4][2]));
   polygon3[4] = G4TwoVector(innercut*sin(ModuleEuler[3][1])*cos(ModuleEuler[3][2]),innercut*sin(ModuleEuler[3][1])*sin(ModuleEuler[3][2]));
   G4ExtrudedSolid* PentaCut = new G4ExtrudedSolid("PentaCut",  polygon3, halfheight, G4TwoVector(0, 0), 0.00001, G4TwoVector(0, 0), 1);
   G4LogicalVolume* logicPentaCut = new G4LogicalVolume(PentaCut, matShell, "Shell_log", 0, 0, 0 );
   G4SubtractionSolid* CutPenta = new G4SubtractionSolid("CutPenta",cutout3, PentaCut, G4Transform3D(NoRotR,G4ThreeVector(0, 0, halfheight)));
   G4SubtractionSolid* CutPenta2 = new G4SubtractionSolid("CutPenta2",cutout30, PentaCut, G4Transform3D(NoRotR,G4ThreeVector(0, 0, halfheight)));
   G4RotationMatrix RotCut = G4RotationMatrix::IDENTITY;
   RotCut.rotateZ(36*degree);
   RotCut.rotateY(180*degree);
   G4SubtractionSolid* DualPentaCut = new G4SubtractionSolid("DualPentaCut", CutPenta, PentaCut, G4Transform3D(RotCut,G4ThreeVector(0, 0, -halfheight)));
  G4SubtractionSolid* DualPentaCut2 = new G4SubtractionSolid("DualPentaCut2", CutPenta2, PentaCut, G4Transform3D(RotCut,G4ThreeVector(0, 0, -halfheight)));
   G4LogicalVolume* logicDualPenta = new G4LogicalVolume(DualPentaCut, matShell, "Shell_log", 0, 0, 0 );
   G4LogicalVolume* logicDualPenta2 = new G4LogicalVolume(DualPentaCut2, matShell, "Shell_log", 0, 0, 0 );
   G4AssemblyVolume* LeftHemi = new G4AssemblyVolume();
   LeftHemi->AddPlacedVolume(logicDualPenta2, NoShiftL, &NoRotL);
   LeftHemi->AddPlacedVolume(logicbump23, NoShiftL, &NoRotL);
   //LeftHemi->AddPlacedVolume(logicbump3, NoShiftL, &NoRotL);
   LeftHemi->MakeImprint(theDetector->HallLog(), NoShiftL, &NoRotL, 0);
   /*
   G4AssemblyVolume* RightHemi = new G4AssemblyVolume();
   RightHemi->AddPlacedVolume(logicDualPenta, NoShiftR, &NoRotR);
   RightHemi->AddPlacedVolume(logicbump10, NoShiftR, &NoRotR);
   //RightHemi->AddPlacedVolume(logicbump30, NoShiftR, &NoRotR);
   RightHemi->MakeImprint(theDetector->HallLog(), NoShiftR, &NoRotR, 0);
   */
}
//Creates and returns the full Greta shell sphere
G4SubtractionSolid* Greta_Shell::Shell(G4String half)
{
  G4double Phi0=0;
  G4double dPhi=0;
  Rot = G4RotationMatrix::IDENTITY;

  if( FindMaterials() ) return NULL;
 
  //////////////////////////////////////////////////
  // The GRETA mounting shell
  //////////////////////////////////////////////////
  if(half == "LEFT"){
    Phi0 = -90.*deg;
    dPhi = 180.*deg;
  }
  else if(half == "RIGHT"){
    Phi0 =  90.*deg;
    dPhi = 180.*deg;
  }
  else if(half == "FULL"){
    Phi0 = 0.*deg;
    dPhi = 360.*deg;
  } else
    G4Exception("Greta_Shell::Shell()", "Error", FatalException, "half argument must be set to LEFT or RIGHT or FULL");

  G4Sphere* solidShell = new G4Sphere( "solidShell", Rmin, Rmax, Phi0, dPhi, 0., 180.*deg);

  // Beam Port  
  G4Tubs* beamPort = new G4Tubs( "beamPort", 0., smallPortRadius, 1.1*Rmax, 0.*deg, 360.*deg);

  G4SubtractionSolid* shell = new G4SubtractionSolid("Shell", solidShell, beamPort, G4Transform3D(Rot,G4ThreeVector(0.0,0.0,0.0)));

  // Subtract the small ports from the shell.
  G4Tubs* smallPort = new G4Tubs("smallPort", 0, smallPortRadius, 1.5*(Rmax-Rmin)/2., 0., 360.*deg);

  G4int iMin = 0;
  G4int iMax = 10;

  for(G4int i = iMin; i<iMax; i++){
    Rot = G4RotationMatrix::IDENTITY;
    Rot.rotateY( PosSP[i].getTheta() );
    Rot.rotateZ( PosSP[i].getPhi() );
    //    if( ( SmallPortStatus[i] <= 0 ) ||( SmallPortStatus[i] >= 0 ) ){
    if( abs(SmallPortStatus[i] <= 1 )) {
      shell = new G4SubtractionSolid ("Shell",shell,smallPort,G4Transform3D(Rot,PosSP[i]));
    }
  }

  // Make space for the LH target
  /*if( half == "LH_north" || half == "LH_south" ){
    G4Box* PortBox = new G4Box("PortBox", modulePortRadius, Rmax-Rmin, modulePortRadius);

    Rot = G4RotationMatrix::IDENTITY;

    shell = new G4SubtractionSolid ("Shell",shell,PortBox,G4Transform3D(Rot,G4ThreeVector(0.,(Rmax+Rmin)/2.,0.0)));
    }*/

  // Subtract the module ports from the shell.
  G4Tubs* modulePort = new G4Tubs("modulePort", 0, modulePortRadius, 1.5*(Rmax-Rmin)/2., 0., 360.*deg);

  iMin = 0;
  iMax = 30;

  for(G4int i = iMin; i<iMax; i++){
    Pos = Pos0;
    Pos.rotateZ(ModuleEuler[i][0]);
    Pos.rotateY(ModuleEuler[i][1]);
    Pos.rotateZ(ModuleEuler[i][2]);
    Rot = G4RotationMatrix::IDENTITY;
    Rot.rotateY( Pos.getTheta() );
    Rot.rotateZ( Pos.getPhi() );

    if( ( ModulePortStatus[i] <= 0 ) || ( ModulePortStatus[i] >= 0 ) ) {
          shell = new G4SubtractionSolid ("Shell",shell,modulePort,G4Transform3D(Rot,Pos));
    }
  }
  return shell;

}
