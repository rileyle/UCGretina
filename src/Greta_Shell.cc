#include "Greta_Shell.hh"

Greta_Shell::Greta_Shell()
{

  // Dimensions taken from fabricationprint_25j2106a-1.pdf
  Rmin      = 1022.0/2.0*mm;
  Rmax      = 1276.0/2.0*mm;

  smallPortRadius  = 7.0/2.0*2.54*cm;   // DETAIL C, SHEET2
  modulePortRadius = 12.008/2.0*2.54*cm; // DETAIL A, SHEET2
  leftOffset       = 0.;
  rightOffset      = 0.;
  G4double RBar = (Rmax+Rmin)/2;
  G4double h = 243.697*mm;
  
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

  // North (Left): -1,  Split: 0,  South (Right): 1, Omit: 2
  SmallPortStatus[0] =  2;
  SmallPortStatus[1] =  2;
  SmallPortStatus[2] =  1;
  SmallPortStatus[3] =  1;
  SmallPortStatus[4] =  0;
  SmallPortStatus[5] =  0;
  SmallPortStatus[6] = -1;
  SmallPortStatus[7] = -1;
  SmallPortStatus[8] =  2;
  SmallPortStatus[9] =  2;

  // For now, assuming hexagons,until Heather Crawford gives us accurate positions.
  //Cutout for Flats  positions
  G4double hPrime = h*(RBar/Rmax);
  G4double d = std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2))*tan(31.717*degree);
  
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

  // North (Left): -1,  Split: 0,  South (Right): 1
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

//Modifies and places the peices of the Greta Shell
void Greta_Shell::Placement(G4String shellStatus,
			    G4bool forwardShellStatus,
			    G4bool backwardShellStatus)
{
  if ( shellStatus != "full" && 
       shellStatus != "left" && 
       shellStatus != "right" ){
    G4cout << "GRETA Shell status " << shellStatus << " is not defined."
	   << G4endl;
    return;
  }

   // Using LTM2_2 - LTM2_5 to find the plane of the flat
  G4double e1 = 0;
  G4double e2 = 31.717473*deg;
  G4double e3 = 18.000000*deg;

  // Geant4 x = Drawing y; Geant4 y = - Drawing x
  G4ThreeVector CADPosFlats[8];
  CADPosFlats[0] = G4ThreeVector(270.000,   87.728, 459.351); // LTM2_0
  CADPosFlats[1] = G4ThreeVector(237.173,  -77.062, 509.936); // LTM2_1
  CADPosFlats[2] = G4ThreeVector(330.754, -107.476, 460.736); // LTM2_2
  CADPosFlats[3] = G4ThreeVector(456.809,   66.984, 353.509); // LTM2_3
  CADPosFlats[4] = G4ThreeVector(389.558,  220.803, 363.662); // LTM2_4
  CADPosFlats[5] = G4ThreeVector(204.406,  281.366, 460.737); // LTM2_5
  CADPosFlats[6] = G4ThreeVector(146.581,  201.751, 509.936); // LTM2_6
  CADPosFlats[7] = G4ThreeVector(0,0,626.3); //center of forward beam port
 
  for(G4int i=0; i<8; i++){
    CADPosFlats[i].rotateZ(-e3);
    CADPosFlats[i].rotateY(-e2);
    CADPosFlats[i].rotateZ(-e1);
  }

  G4double e231 = -90.000000*deg;
  G4double e232 = 121.717473*deg;
  G4double e233 = 90.000000*deg;

  // Geant4 x = Drawing y; Geant4 y = - Drawing x
  G4ThreeVector CADPosHex[4];
  CADPosHex[0] = G4ThreeVector(1.57,   626.5, -120.58); // L23_0
  CADPosHex[1] = G4ThreeVector(-200.94,  560.64, -228.81); // L23_1
  CADPosHex[2] = G4ThreeVector(-219.62, 457.68, -386.44); // L23_2
  CADPosHex[3] = G4ThreeVector(0.41,   388.95, -505.73); // L23_3

  for(G4int i=0; i<4; i++){
    CADPosHex[i].rotateZ(-e233);
    CADPosHex[i].rotateY(-e232);
    CADPosHex[i].rotateZ(-e231);
  }

  G4double e31 = 0*deg;
  G4double e32 = 31.717473*deg;
  G4double e33 = 90.000000*deg;

  // Geant4 x = Drawing y; Geant4 y = - Drawing x
  G4ThreeVector CADPosTriplet[4];
  CADPosTriplet[0] = G4ThreeVector(0,   0, 626.3); // LT_0
  CADPosTriplet[1] = G4ThreeVector(-224.635,  315.287, 507.126); // LT_1
  CADPosTriplet[2] = G4ThreeVector(-101.364, 484.47, 402.565); // LT_2
  CADPosTriplet[3] = G4ThreeVector(87.199, 500.329, 386.149); // LT_3
 
  for(G4int i=0; i<4; i++){
    CADPosTriplet[i].rotateZ(-e33);
    CADPosTriplet[i].rotateY(-e32);
    CADPosTriplet[i].rotateZ(-e31);
  }
  
  G4RunManager* runManager = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector = (DetectorConstruction*) runManager->GetUserDetectorConstruction();
  G4double halfheight = 1300/4*mm;
  G4double RBar = (Rmax+Rmin)/2;
  G4double h = 243.697*mm;
  G4double hPrime = h*(RBar/Rmax);
  G4SubtractionSolid *shellL, *shellR, *shellF;
  std::vector<G4TwoVector> polygon6(5);
  polygon6[0] = G4TwoVector(CADPosFlats[2].x(),CADPosFlats[2].y());
  polygon6[1] = G4TwoVector(CADPosFlats[3].x(),CADPosFlats[3].y());
  polygon6[2] = G4TwoVector(CADPosFlats[4].x(),CADPosFlats[4].y());
  polygon6[3] = G4TwoVector(CADPosFlats[5].x(),CADPosFlats[5].y());
  polygon6[4] = G4TwoVector(CADPosFlats[7].x(),CADPosFlats[7].y());
  std::vector<G4TwoVector> polygon1(4);
  polygon1[0] = G4TwoVector(CADPosHex[0].x(),CADPosHex[0].y());
  polygon1[1] = G4TwoVector(CADPosHex[1].x(),CADPosHex[1].y());
  polygon1[2] = G4TwoVector(CADPosHex[2].x(),CADPosHex[2].y());
  polygon1[3] = G4TwoVector(CADPosHex[3].x(),CADPosHex[3].y());
  std::vector<G4TwoVector> polygon2(4);
  polygon2[0] = G4TwoVector(CADPosTriplet[0].x(),CADPosTriplet[0].y());
  polygon2[1] = G4TwoVector(CADPosTriplet[1].x(),CADPosTriplet[1].y());
  polygon2[2] = G4TwoVector(CADPosTriplet[2].x(),CADPosTriplet[2].y());
  polygon2[3] = G4TwoVector(CADPosTriplet[3].x(),CADPosTriplet[3].y());
  
  std::vector<G4ExtrudedSolid::ZSection> zsections;
  zsections.push_back(G4ExtrudedSolid::ZSection(0,G4TwoVector(0,0), 1.1));
  zsections.push_back(G4ExtrudedSolid::ZSection((Rmin+Rmax)/2,G4TwoVector(0,0), 3));
  G4ExtrudedSolid* tripleShape = new G4ExtrudedSolid("tripleShape",  polygon6, zsections);
  
  if(shellStatus != "full"){
    shellL = Shell("LEFT");
    shellR = Shell("RIGHT");
  }
  else if(shellStatus == "full"){
    shellF = Shell("FULL");
  }
 
		   
  G4double halfside = 230*mm;
  G4double innercut = 751*mm;
  G4ThreeVector NoShiftR = G4ThreeVector(0, 0, 0);
  G4RotationMatrix NoRotR = G4RotationMatrix::IDENTITY;
   
  G4ExtrudedSolid* solidTarget1 = new G4ExtrudedSolid("solidTarget1",  polygon1, halfheight, G4TwoVector(0, 0), 0.00001, G4TwoVector(0, 0), 1);

  G4ExtrudedSolid* solidTarget2 = new G4ExtrudedSolid("solidTarget2",  polygon2, halfheight, G4TwoVector(0, 0), 0.00001, G4TwoVector(0, 0), 1);
 
  //HexHole for Hole 10
  G4ThreeVector TarPos10 = G4ThreeVector(0, 0, halfheight);
  TarPos10.rotateZ(ModuleEuler[9][0]);
  TarPos10.rotateY(ModuleEuler[9][1]);
  TarPos10.rotateZ(ModuleEuler[9][2]);
  G4RotationMatrix RotShell10 = G4RotationMatrix::IDENTITY;
  RotShell10.rotateZ( -90*degree );
  RotShell10.rotateY( TarPos10.getTheta() );
  RotShell10.rotateZ( TarPos10.getPhi() );
  G4SubtractionSolid* cutout10 = new G4SubtractionSolid("cutout10", shellL, solidTarget1, G4Transform3D(RotShell10, TarPos10));

  // bumpout for the Right "hemisphere"
  G4IntersectionSolid* bump10 = new G4IntersectionSolid("bump10",shellL, solidTarget1, G4Transform3D(RotShell10, TarPos10));
  G4LogicalVolume* logicbump10 = new G4LogicalVolume(bump10, matShell, "Shell_log", 0, 0, 0 );

  //TripleHole for Hole 30
  G4ThreeVector TarPos30 = G4ThreeVector(0, 0, halfheight);
  TarPos30.rotateZ(ModuleEuler[29][0]);
  TarPos30.rotateY(ModuleEuler[29][1]);
  TarPos30.rotateZ(ModuleEuler[29][2]);
  G4RotationMatrix RotShell30 = G4RotationMatrix::IDENTITY;
  RotShell30.rotateZ( 7*degree );
  RotShell30.rotateY( TarPos30.getTheta() );
  RotShell30.rotateZ( TarPos30.getPhi() );
  // Left hemisphere with both cutouts but not the bumpouts
  G4SubtractionSolid* cutout30 = new G4SubtractionSolid("cutout30", cutout10, solidTarget2, G4Transform3D(RotShell30, TarPos30));
  G4LogicalVolume* logicCutout30 = new G4LogicalVolume(cutout30, matShell, "Shell_log", 0, 0, 0 );
     
  // bumpout for the Right "hemisphere"
  G4IntersectionSolid* bump30 = new G4IntersectionSolid("bump30",cutout10, solidTarget2, G4Transform3D(RotShell30, TarPos30));
  G4LogicalVolume* logicbump30 = new G4LogicalVolume(bump30, matShell, "Shell_log", 0, 0, 0 );
  
  //HexHole for Hole 23
  G4ThreeVector TarPos23 = G4ThreeVector(0, 0, halfheight);
  TarPos23.rotateZ(ModuleEuler[22][0]);
  TarPos23.rotateY(ModuleEuler[22][1]);
  TarPos23.rotateZ(ModuleEuler[22][2]);
  G4RotationMatrix RotShell23 = G4RotationMatrix::IDENTITY;
  RotShell23.rotateZ( -90*degree );
  RotShell23.rotateY( TarPos23.getTheta() );
  RotShell23.rotateZ( TarPos23.getPhi() );
  // Right hemisphere with the hole 23 cutout
  G4SubtractionSolid* cutout23 = new G4SubtractionSolid("cutout23", shellR, solidTarget1, G4Transform3D(RotShell23, TarPos23));

  // bumpout for the Left hemisphere.
  G4IntersectionSolid* bump23 = new G4IntersectionSolid("bump23", shellR, solidTarget1, G4Transform3D(RotShell23, TarPos23));
  G4LogicalVolume* logicbump23 = new G4LogicalVolume(bump23, matShell, "Shell_log", 0, 0, 0 );
    
  //TripleHole for Hole 3
  G4ThreeVector TarPos3 = G4ThreeVector(0, 0, halfheight);
  TarPos3.rotateZ(ModuleEuler[2][0]);
  TarPos3.rotateY(ModuleEuler[2][1]);
  TarPos3.rotateZ(ModuleEuler[2][2]);
  G4RotationMatrix RotShell3 = G4RotationMatrix::IDENTITY;
  RotShell3.rotateZ( 7*degree );
  RotShell3.rotateY( TarPos3.getTheta() );
  RotShell3.rotateZ( TarPos3.getPhi() );
  // Right hemisphere with both cutouts but not the bumpouts
  G4SubtractionSolid* cutout3 = new G4SubtractionSolid("cutout3", cutout23, solidTarget2, G4Transform3D(RotShell3, TarPos3));
  G4LogicalVolume* logicCutout3 = new G4LogicalVolume(cutout3, matShell, "Shell_log", 0, 0, 0 );
  // bumpout for the Left hemisphere
  G4IntersectionSolid* bump3 = new G4IntersectionSolid("bump3", cutout23, solidTarget2, G4Transform3D(RotShell3, TarPos3));
  G4LogicalVolume* logicbump3 = new G4LogicalVolume(bump3, matShell, "Shell_log", 0, 0, 0 );

  G4ThreeVector NoShiftL = G4ThreeVector(0, 0, 0);
  G4RotationMatrix NoRotL = G4RotationMatrix::IDENTITY;

  if(!forwardShellStatus || !backwardShellStatus) {
    //ForwardPentaCut polygon
    std::vector<G4TwoVector> polygon3(5);
    polygon3[0] = G4TwoVector(innercut*sin(ModuleEuler[2][1])*cos(ModuleEuler[2][2]),
			      innercut*sin(ModuleEuler[2][1])*sin(ModuleEuler[2][2]));
    polygon3[1] = G4TwoVector(innercut*sin(ModuleEuler[1][1])*cos(ModuleEuler[1][2]),
			      innercut*sin(ModuleEuler[1][1])*sin(ModuleEuler[1][2]));
    polygon3[2] = G4TwoVector(innercut*sin(ModuleEuler[0][1])*cos(ModuleEuler[0][2]),
			      innercut*sin(ModuleEuler[0][1])*sin(ModuleEuler[0][2]));
    polygon3[3] = G4TwoVector(innercut*sin(ModuleEuler[4][1])*cos(ModuleEuler[4][2]),
			      innercut*sin(ModuleEuler[4][1])*sin(ModuleEuler[4][2]));
    polygon3[4] = G4TwoVector(innercut*sin(ModuleEuler[3][1])*cos(ModuleEuler[3][2]),
			      innercut*sin(ModuleEuler[3][1])*sin(ModuleEuler[3][2]));
    G4ExtrudedSolid* PentaCut = new G4ExtrudedSolid("PentaCut",  polygon3, halfheight,
						    G4TwoVector(0, 0), 0.00001, G4TwoVector(0, 0), 1);
    //BackwardPentaCut polygon
    std::vector<G4TwoVector> polygon7(5);
    polygon7[0] = G4TwoVector(innercut*sin(ModuleEuler[29][1])*cos(ModuleEuler[29][2]),
			      innercut*sin(ModuleEuler[29][1])*sin(ModuleEuler[29][2]));
    polygon7[1] = G4TwoVector(innercut*sin(ModuleEuler[28][1])*cos(ModuleEuler[28][2]),
			      innercut*sin(ModuleEuler[28][1])*sin(ModuleEuler[28][2]));
    polygon7[2] = G4TwoVector(innercut*sin(ModuleEuler[27][1])*cos(ModuleEuler[27][2]),
			      innercut*sin(ModuleEuler[27][1])*sin(ModuleEuler[27][2]));
    polygon7[3] = G4TwoVector(innercut*sin(ModuleEuler[26][1])*cos(ModuleEuler[26][2]),
			      innercut*sin(ModuleEuler[26][1])*sin(ModuleEuler[26][2]));
    polygon7[4] = G4TwoVector(innercut*sin(ModuleEuler[25][1])*cos(ModuleEuler[25][2]),
			      innercut*sin(ModuleEuler[25][1])*sin(ModuleEuler[25][2]));
    G4ExtrudedSolid* BPentaCut = new G4ExtrudedSolid("BPentaCut",  polygon7, halfheight,
						     G4TwoVector(0, 0), 1, G4TwoVector(0, 0), 0.00001);
    G4RotationMatrix RotCut = G4RotationMatrix::IDENTITY;
    RotCut.rotateZ(36*degree);
    RotCut.rotateY(180*degree);
    if(shellStatus == "left" && !forwardShellStatus && backwardShellStatus) {
      G4SubtractionSolid* CutPentaLFo
	= new G4SubtractionSolid("CutPentaLFo",cutout30, PentaCut,
				 G4Transform3D(NoRotR, G4ThreeVector(0, 0, halfheight)));

      // Cut the flats for the (omitted) forward triple
      G4ThreeVector TarPos1 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos1.rotateY(ModuleEuler[0][1]);
      TarPos1.rotateZ(ModuleEuler[0][2]);
      G4RotationMatrix RotShell1 = G4RotationMatrix::IDENTITY;
      RotShell1.rotateZ(ModuleEuler[0][0]);
      RotShell1.rotateY( TarPos1.getTheta() );
      RotShell1.rotateZ( TarPos1.getPhi() );

      G4ThreeVector TarPos2 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos2.rotateY(ModuleEuler[1][1]);
      TarPos2.rotateZ(ModuleEuler[1][2]);
      G4RotationMatrix RotShell2 = G4RotationMatrix::IDENTITY;
      RotShell2.rotateZ(ModuleEuler[1][0]);
      RotShell2.rotateY( TarPos2.getTheta() );
      RotShell2.rotateZ( TarPos2.getPhi() );
      
      G4ThreeVector TarPos3 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos3.rotateY(ModuleEuler[2][1]);
      TarPos3.rotateZ(ModuleEuler[2][2]);
      G4RotationMatrix RotShell3 = G4RotationMatrix::IDENTITY;
      RotShell3.rotateZ(ModuleEuler[2][0]);
      RotShell3.rotateY( TarPos3.getTheta() );
      RotShell3.rotateZ( TarPos3.getPhi() );

      G4SubtractionSolid* cutLTriple = new G4SubtractionSolid("cutLTriple", CutPentaLFo, tripleShape, G4Transform3D(RotShell1, TarPos1));
      G4SubtractionSolid* cutLTriple2 = new G4SubtractionSolid("cutLTriple", cutLTriple, tripleShape, G4Transform3D(RotShell2, TarPos2));
      G4SubtractionSolid* cutLTriple3 = new G4SubtractionSolid("cutLTriple", cutLTriple2, tripleShape, G4Transform3D(RotShell3, TarPos3));
      G4LogicalVolume* logiccutLTriple3 = new G4LogicalVolume(cutLTriple3, matShell, "Shell_log", 0, 0, 0 );
      
      //Place Left Hemisphere with full cutouts, bumpouts, and no forward ring
      G4AssemblyVolume* LeftHemi = new G4AssemblyVolume();
      LeftHemi->AddPlacedVolume(logiccutLTriple3, NoShiftR, &NoRotR);
      LeftHemi->AddPlacedVolume(logicbump23, NoShiftR, &NoRotR);
      LeftHemi->MakeImprint(theDetector->HallLog(), NoShiftR, &NoRotR, 0);
    }
    else if(shellStatus == "left" && !backwardShellStatus && forwardShellStatus){
      G4SubtractionSolid* CutPentaLBo
	= new G4SubtractionSolid("CutPentaLBo",cutout30, PentaCut,
				 G4Transform3D(RotCut, G4ThreeVector(0, 0, -halfheight)));
      G4ThreeVector TarPos26 = G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos26.rotateY( ModuleEuler[25][1] );
      TarPos26.rotateZ( ModuleEuler[25][2] );
      G4RotationMatrix RotShell26 = G4RotationMatrix::IDENTITY;
      RotShell26.rotateZ( ModuleEuler[25][0] ); // Orient the polygon in the plane tangent to the sphere
      RotShell26.rotateY( TarPos26.getTheta() );
      RotShell26.rotateZ( TarPos26.getPhi() );
      G4SubtractionSolid* LDoubleFlat = new G4SubtractionSolid("LDoubleFlat", CutPentaLBo, tripleShape, G4Transform3D(RotShell26, TarPos26));
      G4ThreeVector TarPos27 = G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos27.rotateY( ModuleEuler[26][1] );
      TarPos27.rotateZ( ModuleEuler[26][2] );
      G4RotationMatrix RotShell27 = G4RotationMatrix::IDENTITY;
      RotShell27.rotateZ( ModuleEuler[26][0] ); // Orient the polygon in the plane tangent to the sphere
      RotShell27.rotateY( TarPos27.getTheta() );
      RotShell27.rotateZ( TarPos27.getPhi() );
      G4SubtractionSolid* LDoubleFlat2 = new G4SubtractionSolid("LDoubleFlat2", LDoubleFlat, tripleShape, G4Transform3D(RotShell27, TarPos27));
      G4LogicalVolume* logicLDoubleFlat2 = new G4LogicalVolume(LDoubleFlat2, matShell, "Shell_log", 0, 0, 0 );
      
      //Place Left Hemisphere with full cutouts, bumpouts, and no backward ring
      G4AssemblyVolume* LeftHemi = new G4AssemblyVolume();
      LeftHemi->AddPlacedVolume(logicLDoubleFlat2, NoShiftR, &NoRotR);
      LeftHemi->AddPlacedVolume(logicbump23, NoShiftR, &NoRotR);
      LeftHemi->MakeImprint(theDetector->HallLog(), NoShiftR, &NoRotR, 0);
    }
    else if(shellStatus == "left" && !backwardShellStatus && !forwardShellStatus) {
      G4SubtractionSolid* CutPentaLFB
	= new G4SubtractionSolid("CutPentaLF", cutout30, PentaCut,
				 G4Transform3D(NoRotR, G4ThreeVector(0,0,halfheight)));
      G4SubtractionSolid* DualPentaCutLFB
	= new G4SubtractionSolid("DualPentaCutLFB", CutPentaLFB, PentaCut,
				 G4Transform3D(RotCut, G4ThreeVector(0, 0, -halfheight)));
      G4ThreeVector TarPos26 = G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos26.rotateY( ModuleEuler[25][1] );
      TarPos26.rotateZ( ModuleEuler[25][2] );
      G4RotationMatrix RotShell26 = G4RotationMatrix::IDENTITY;
      RotShell26.rotateZ( ModuleEuler[25][0] ); // Orient the polygon in the plane tangent to the sphere
      RotShell26.rotateY( TarPos26.getTheta() );
      RotShell26.rotateZ( TarPos26.getPhi() );
      G4SubtractionSolid* LDoubleFlat = new G4SubtractionSolid("LDoubleFlat", DualPentaCutLFB, tripleShape, G4Transform3D(RotShell26, TarPos26));
      G4ThreeVector TarPos27 = G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos27.rotateY( ModuleEuler[26][1] );
      TarPos27.rotateZ( ModuleEuler[26][2] );
      G4RotationMatrix RotShell27 = G4RotationMatrix::IDENTITY;
      RotShell27.rotateZ( ModuleEuler[26][0] ); // Orient the polygon in the plane tangent to the sphere
      RotShell27.rotateY( TarPos27.getTheta() );
      RotShell27.rotateZ( TarPos27.getPhi() );
      G4SubtractionSolid* LDoubleFlat2 = new G4SubtractionSolid("LDoubleFlat2", LDoubleFlat, tripleShape, G4Transform3D(RotShell27, TarPos27));

      G4ThreeVector TarPos1 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos1.rotateY(ModuleEuler[0][1]);
      TarPos1.rotateZ(ModuleEuler[0][2]);
      G4RotationMatrix RotShell1 = G4RotationMatrix::IDENTITY;
      RotShell1.rotateZ(ModuleEuler[0][0]);
      RotShell1.rotateY( TarPos1.getTheta() );
      RotShell1.rotateZ( TarPos1.getPhi() );

      G4ThreeVector TarPos2 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos2.rotateY(ModuleEuler[1][1]);
      TarPos2.rotateZ(ModuleEuler[1][2]);
      G4RotationMatrix RotShell2 = G4RotationMatrix::IDENTITY;
      RotShell2.rotateZ(ModuleEuler[1][0]);
      RotShell2.rotateY( TarPos2.getTheta() );
      RotShell2.rotateZ( TarPos2.getPhi() );
      
      G4ThreeVector TarPos3 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos3.rotateY(ModuleEuler[2][1]);
      TarPos3.rotateZ(ModuleEuler[2][2]);
      G4RotationMatrix RotShell3 = G4RotationMatrix::IDENTITY;
      RotShell3.rotateZ(ModuleEuler[2][0]);
      RotShell3.rotateY( TarPos3.getTheta() );
      RotShell3.rotateZ( TarPos3.getPhi() );

      G4SubtractionSolid* LTripleFlat = new G4SubtractionSolid("LTripleFlat", LDoubleFlat2, tripleShape, G4Transform3D(RotShell1, TarPos1));
      G4SubtractionSolid* LTripleFlat2 = new G4SubtractionSolid("LTripleFlat2", LTripleFlat, tripleShape, G4Transform3D(RotShell2, TarPos2));
      G4SubtractionSolid* LTripleFlat3 = new G4SubtractionSolid("LTripleFlat3", LTripleFlat2, tripleShape, G4Transform3D(RotShell3, TarPos3));
      G4LogicalVolume* logicLTripleFlat3 = new G4LogicalVolume(LTripleFlat3, matShell, "Shell_log", 0, 0, 0 );
      
      //Place Left Hemisphere with full cutouts, bumpouts and no forward or backward ring
      G4AssemblyVolume* LeftHemi = new G4AssemblyVolume();
      LeftHemi->AddPlacedVolume(logicLTripleFlat3, NoShiftL, &NoRotL);
      LeftHemi->AddPlacedVolume(logicbump23, NoShiftL, &NoRotL);
      LeftHemi->MakeImprint(theDetector->HallLog(), NoShiftL, &NoRotL, 0);
    }
    if(shellStatus == "right" && !forwardShellStatus && backwardShellStatus) {
      G4SubtractionSolid* CutPentaRFo
	= new G4SubtractionSolid("CutPentaRFo",cutout3, PentaCut,
				 G4Transform3D(NoRotR, G4ThreeVector(0, 0, halfheight)));
      G4ThreeVector TarPos4 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos4.rotateY(ModuleEuler[3][1]);
      TarPos4.rotateZ(ModuleEuler[3][2]);
      G4RotationMatrix RotShell4 = G4RotationMatrix::IDENTITY;
      RotShell4.rotateZ(ModuleEuler[3][0]);
      RotShell4.rotateY( TarPos4.getTheta() );
      RotShell4.rotateZ( TarPos4.getPhi() );
      G4SubtractionSolid* RDoubleFlat = new G4SubtractionSolid("RDoubleFlat", CutPentaRFo, tripleShape, G4Transform3D(RotShell4, TarPos4));
      G4ThreeVector TarPos5 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos5.rotateY(ModuleEuler[4][1]);
      TarPos5.rotateZ(ModuleEuler[4][2]);
      G4RotationMatrix RotShell5 = G4RotationMatrix::IDENTITY;
      RotShell5.rotateZ(ModuleEuler[4][0]);
      RotShell5.rotateY( TarPos5.getTheta() );
      RotShell5.rotateZ( TarPos5.getPhi() );
      G4SubtractionSolid* RDoubleFlat2 = new G4SubtractionSolid("RDoubleFlat2", RDoubleFlat, tripleShape, G4Transform3D(RotShell5, TarPos5));
      G4LogicalVolume* logicRDoubleFlat2 = new G4LogicalVolume(RDoubleFlat2, matShell, "Shell_log", 0, 0, 0 );
      
      //Place Right Hemisphere with full cutouts, bumpouts, and no forward ring
      G4AssemblyVolume* RightHemi = new G4AssemblyVolume();
      RightHemi->AddPlacedVolume(logicRDoubleFlat2, NoShiftR, &NoRotR);
      RightHemi->AddPlacedVolume(logicbump10, NoShiftR, &NoRotR);
      RightHemi->MakeImprint(theDetector->HallLog(), NoShiftR, &NoRotR, 0);
    }
    else if(shellStatus == "right" && !backwardShellStatus && forwardShellStatus) {
      G4SubtractionSolid* CutPentaRBo
	= new G4SubtractionSolid("CutPentaRBo",cutout3, PentaCut,
				 G4Transform3D(RotCut, G4ThreeVector(0, 0, -halfheight)));
      G4ThreeVector TarPos28 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos28.rotateY(ModuleEuler[27][1]);
      TarPos28.rotateZ(ModuleEuler[27][2]);
      G4RotationMatrix RotShell28 = G4RotationMatrix::IDENTITY;
      RotShell28.rotateZ(ModuleEuler[27][0]);
      RotShell28.rotateY( TarPos28.getTheta() );
      RotShell28.rotateZ( TarPos28.getPhi() );
      G4SubtractionSolid* RTripleFlat = new G4SubtractionSolid("RTripleFlat", CutPentaRBo, tripleShape, G4Transform3D(RotShell28, TarPos28));
      G4ThreeVector TarPos29 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos29.rotateY(ModuleEuler[28][1]);
      TarPos29.rotateZ(ModuleEuler[28][2]);
      G4RotationMatrix RotShell29 = G4RotationMatrix::IDENTITY;
      RotShell29.rotateZ(ModuleEuler[28][0]);
      RotShell29.rotateY( TarPos29.getTheta() );
      RotShell29.rotateZ( TarPos29.getPhi() );
      G4SubtractionSolid* RTripleFlat2 = new G4SubtractionSolid("RTripleFlat2", RTripleFlat, tripleShape, G4Transform3D(RotShell29, TarPos29));
      G4ThreeVector TarPos30 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos30.rotateY(ModuleEuler[29][1]);
      TarPos30.rotateZ(ModuleEuler[29][2]);
      G4RotationMatrix RotShell30 = G4RotationMatrix::IDENTITY;
      RotShell30.rotateZ(ModuleEuler[29][0]);
      RotShell30.rotateY( TarPos30.getTheta() );
      RotShell30.rotateZ( TarPos30.getPhi() );
      G4SubtractionSolid* RTripleFlat3 = new G4SubtractionSolid("RTripleFlat3", RTripleFlat2, tripleShape, G4Transform3D(RotShell30, TarPos30));
      G4LogicalVolume* logicRTripleFlat3 = new G4LogicalVolume(RTripleFlat3, matShell, "Shell_log", 0, 0, 0 );
      
      //Place Right Hemisphere with full cutouts, bumpouts, and no backward ring
      G4AssemblyVolume* RightHemi = new G4AssemblyVolume();
      RightHemi->AddPlacedVolume(logicRTripleFlat3, NoShiftR, &NoRotR);
      RightHemi->AddPlacedVolume(logicbump10, NoShiftR, &NoRotR);
      RightHemi->MakeImprint(theDetector->HallLog(), NoShiftR, &NoRotR, 0);
    }
    else if(shellStatus == "right" && !backwardShellStatus && !forwardShellStatus) {
      G4SubtractionSolid* CutPentaRF
	= new G4SubtractionSolid("CutPentaRF", cutout3, PentaCut,
				 G4Transform3D(NoRotR, G4ThreeVector(0, 0, halfheight)));
      G4SubtractionSolid* DualPentaCutRFB
	= new G4SubtractionSolid("DualPentaCutRFB", CutPentaRF, PentaCut,
				 G4Transform3D(RotCut, G4ThreeVector(0, 0, -halfheight)));
      
      G4ThreeVector TarPos4 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos4.rotateY(ModuleEuler[3][1]);
      TarPos4.rotateZ(ModuleEuler[3][2]);
      G4RotationMatrix RotShell4 = G4RotationMatrix::IDENTITY;
      RotShell4.rotateZ(ModuleEuler[3][0]);
      RotShell4.rotateY( TarPos4.getTheta() );
      RotShell4.rotateZ( TarPos4.getPhi() );
      G4SubtractionSolid* RDoubleFlat = new G4SubtractionSolid("RDoubleFlat", DualPentaCutRFB, tripleShape, G4Transform3D(RotShell4, TarPos4));
      G4ThreeVector TarPos5 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos5.rotateY(ModuleEuler[4][1]);
      TarPos5.rotateZ(ModuleEuler[4][2]);
      G4RotationMatrix RotShell5 = G4RotationMatrix::IDENTITY;
      RotShell5.rotateZ(ModuleEuler[4][0]);
      RotShell5.rotateY( TarPos5.getTheta() );
      RotShell5.rotateZ( TarPos5.getPhi() );
      G4SubtractionSolid* RDoubleFlat2 = new G4SubtractionSolid("RDoubleFlat2", RDoubleFlat, tripleShape, G4Transform3D(RotShell5, TarPos5));

      G4ThreeVector TarPos28 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos28.rotateY(ModuleEuler[27][1]);
      TarPos28.rotateZ(ModuleEuler[27][2]);
      G4RotationMatrix RotShell28 = G4RotationMatrix::IDENTITY;
      RotShell28.rotateZ(ModuleEuler[27][0]);
      RotShell28.rotateY( TarPos28.getTheta() );
      RotShell28.rotateZ( TarPos28.getPhi() );
      G4SubtractionSolid* RTripleFlat = new G4SubtractionSolid("RTripleFlat", RDoubleFlat2, tripleShape, G4Transform3D(RotShell28, TarPos28));
      G4ThreeVector TarPos29 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos29.rotateY(ModuleEuler[28][1]);
      TarPos29.rotateZ(ModuleEuler[28][2]);
      G4RotationMatrix RotShell29 = G4RotationMatrix::IDENTITY;
      RotShell29.rotateZ(ModuleEuler[28][0]);
      RotShell29.rotateY( TarPos29.getTheta() );
      RotShell29.rotateZ( TarPos29.getPhi() );
      G4SubtractionSolid* RTripleFlat2 = new G4SubtractionSolid("RTripleFlat2", RTripleFlat, tripleShape, G4Transform3D(RotShell29, TarPos29));
      G4ThreeVector TarPos30 =  G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos30.rotateY(ModuleEuler[29][1]);
      TarPos30.rotateZ(ModuleEuler[29][2]);
      G4RotationMatrix RotShell30 = G4RotationMatrix::IDENTITY;
      RotShell30.rotateZ(ModuleEuler[29][0]);
      RotShell30.rotateY( TarPos30.getTheta() );
      RotShell30.rotateZ( TarPos30.getPhi() );
      G4SubtractionSolid* RTripleFlat3 = new G4SubtractionSolid("RTripleFlat3", RTripleFlat2, tripleShape, G4Transform3D(RotShell30, TarPos30));
      G4LogicalVolume* logicRTripleFlat3 = new G4LogicalVolume(RTripleFlat3, matShell, "Shell_log", 0, 0, 0 );
      
      //Place Right Hemisphere with full cutouts, bumpouts, and no forward or backward ring
      G4AssemblyVolume* RightHemi = new G4AssemblyVolume();
      RightHemi->AddPlacedVolume(logicRTripleFlat3, NoShiftR, &NoRotR);
      RightHemi->AddPlacedVolume(logicbump10, NoShiftR, &NoRotR);
      RightHemi->MakeImprint(theDetector->HallLog(), NoShiftR, &NoRotR, 0);
    }
    if(shellStatus == "full" && !forwardShellStatus && backwardShellStatus) {

      // Define the transformations for the flats.
      // (Eventually, this should be a for loop.)
      G4ThreeVector TarPos1
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos1.rotateY(ModuleEuler[0][1]);
      TarPos1.rotateZ(ModuleEuler[0][2]);
      G4RotationMatrix RotShell1 = G4RotationMatrix::IDENTITY;
      RotShell1.rotateZ(ModuleEuler[0][0]);
      RotShell1.rotateY( TarPos1.getTheta() );
      RotShell1.rotateZ( TarPos1.getPhi() );

      G4ThreeVector TarPos2
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos2.rotateY(ModuleEuler[1][1]);
      TarPos2.rotateZ(ModuleEuler[1][2]);
      G4RotationMatrix RotShell2 = G4RotationMatrix::IDENTITY;
      RotShell2.rotateZ(ModuleEuler[1][0]);
      RotShell2.rotateY( TarPos2.getTheta() );
      RotShell2.rotateZ( TarPos2.getPhi() );
      
      G4ThreeVector TarPos3
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos3.rotateY(ModuleEuler[2][1]);
      TarPos3.rotateZ(ModuleEuler[2][2]);
      G4RotationMatrix RotShell3 = G4RotationMatrix::IDENTITY;
      RotShell3.rotateZ(ModuleEuler[2][0]);
      RotShell3.rotateY( TarPos3.getTheta() );
      RotShell3.rotateZ( TarPos3.getPhi() );

      G4ThreeVector TarPos4
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos4.rotateY(ModuleEuler[3][1]);
      TarPos4.rotateZ(ModuleEuler[3][2]);
      G4RotationMatrix RotShell4 = G4RotationMatrix::IDENTITY;
      RotShell4.rotateZ(ModuleEuler[3][0]);
      RotShell4.rotateY( TarPos4.getTheta() );
      RotShell4.rotateZ( TarPos4.getPhi() );

      G4ThreeVector TarPos5
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos5.rotateY(ModuleEuler[4][1]);
      TarPos5.rotateZ(ModuleEuler[4][2]);
      G4RotationMatrix RotShell5 = G4RotationMatrix::IDENTITY;
      RotShell5.rotateZ(ModuleEuler[4][0]);
      RotShell5.rotateY( TarPos5.getTheta() );
      RotShell5.rotateZ( TarPos5.getPhi() );

      // Cut the flats first ...
      G4SubtractionSolid* FTripleFlat
	= new G4SubtractionSolid("FTripleFlat", shellF, tripleShape,
				 G4Transform3D(RotShell1, TarPos1));
      G4SubtractionSolid* FTripleFlat2
	= new G4SubtractionSolid("FTripleFlat2", FTripleFlat, tripleShape,
				 G4Transform3D(RotShell2, TarPos2));
      G4SubtractionSolid* FTripleFlat3
	= new G4SubtractionSolid("FTripleFlat3", FTripleFlat2, tripleShape,
				 G4Transform3D(RotShell3, TarPos3));
      G4SubtractionSolid* FDoubleFlat
	= new G4SubtractionSolid("FDoubleFlat", FTripleFlat3, tripleShape,
				 G4Transform3D(RotShell4, TarPos4));
      G4SubtractionSolid* FDoubleFlat2
	= new G4SubtractionSolid("FDoubleFlat2", FDoubleFlat, tripleShape,
				 G4Transform3D(RotShell5, TarPos5));

      // ... then remove the central pentagon.
      G4SubtractionSolid* FFRingCutPenta
	= new G4SubtractionSolid("FRingFlat", FDoubleFlat2, PentaCut,
				 G4Transform3D(NoRotR,
					       G4ThreeVector(0,0,halfheight)));

      G4LogicalVolume* logicFFRingCutPenta
	= new G4LogicalVolume(FFRingCutPenta, matShell, "Shell_log", 0, 0, 0 );
      
      //Place Full Shell with no forward ring
      G4AssemblyVolume* FullNF = new G4AssemblyVolume();
      FullNF->AddPlacedVolume(logicFFRingCutPenta, NoShiftR, &NoRotR);
      FullNF->MakeImprint(theDetector->HallLog(), NoShiftR, &NoRotR, 0);
    }
    else if(shellStatus == "full" && !backwardShellStatus && forwardShellStatus) {

      // Define the transformations for the flats.
      // (Eventually, this should be a for loop.)
      G4ThreeVector TarPos26
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos26.rotateY(ModuleEuler[25][1]);
      TarPos26.rotateZ(ModuleEuler[25][2]);
      G4RotationMatrix RotShell26 = G4RotationMatrix::IDENTITY;
      RotShell26.rotateZ(ModuleEuler[25][0]);
      RotShell26.rotateY( TarPos26.getTheta() );
      RotShell26.rotateZ( TarPos26.getPhi() );

      G4ThreeVector TarPos27
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos27.rotateY(ModuleEuler[26][1]);
      TarPos27.rotateZ(ModuleEuler[26][2]);
      G4RotationMatrix RotShell27 = G4RotationMatrix::IDENTITY;
      RotShell27.rotateZ(ModuleEuler[26][0]);
      RotShell27.rotateY( TarPos27.getTheta() );
      RotShell27.rotateZ( TarPos27.getPhi() );
      
      G4ThreeVector TarPos28
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos28.rotateY(ModuleEuler[27][1]);
      TarPos28.rotateZ(ModuleEuler[27][2]);
      G4RotationMatrix RotShell28 = G4RotationMatrix::IDENTITY;
      RotShell28.rotateZ(ModuleEuler[27][0]);
      RotShell28.rotateY( TarPos28.getTheta() );
      RotShell28.rotateZ( TarPos28.getPhi() );

      G4ThreeVector TarPos29
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos29.rotateY(ModuleEuler[28][1]);
      TarPos29.rotateZ(ModuleEuler[28][2]);
      G4RotationMatrix RotShell29 = G4RotationMatrix::IDENTITY;
      RotShell29.rotateZ(ModuleEuler[28][0]);
      RotShell29.rotateY( TarPos29.getTheta() );
      RotShell29.rotateZ( TarPos29.getPhi() );

      G4ThreeVector TarPos30
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos30.rotateY(ModuleEuler[29][1]);
      TarPos30.rotateZ(ModuleEuler[29][2]);
      G4RotationMatrix RotShell30 = G4RotationMatrix::IDENTITY;
      RotShell30.rotateZ(ModuleEuler[29][0]);
      RotShell30.rotateY( TarPos30.getTheta() );
      RotShell30.rotateZ( TarPos30.getPhi() );

      // Cut the flats first ...
      G4SubtractionSolid* BTripleFlat
	= new G4SubtractionSolid("BTripleFlat", shellF, tripleShape,
				 G4Transform3D(RotShell26, TarPos26));
      G4SubtractionSolid* BTripleFlat2
	= new G4SubtractionSolid("BTripleFlat2", BTripleFlat, tripleShape,
				 G4Transform3D(RotShell27, TarPos27));
      G4SubtractionSolid* BTripleFlat3
	= new G4SubtractionSolid("BTripleFlat3", BTripleFlat2, tripleShape,
				 G4Transform3D(RotShell28, TarPos28));
      G4SubtractionSolid* BDoubleFlat
	= new G4SubtractionSolid("BDoubleFlat", BTripleFlat3, tripleShape,
				 G4Transform3D(RotShell29, TarPos29));
      G4SubtractionSolid* BDoubleFlat2
	= new G4SubtractionSolid("BDoubleFlat2", BDoubleFlat, tripleShape,
				 G4Transform3D(RotShell30, TarPos30));

      // ... then remove the central pentagon (rotated for the back ring).
      G4SubtractionSolid* FBRingCutPenta
	= new G4SubtractionSolid("BRingCutPenta", BDoubleFlat2, BPentaCut,
				 G4Transform3D(NoRotR,
					       G4ThreeVector(0,0,-halfheight)));

      G4LogicalVolume* logicFBRingCutPenta
	= new G4LogicalVolume(FBRingCutPenta, matShell, "Shell_log", 0, 0, 0 );
      
      //Place Full Shell with no backward ring
      G4AssemblyVolume* FullNB = new G4AssemblyVolume();
      FullNB->AddPlacedVolume(logicFBRingCutPenta, NoShiftR, &NoRotR);
      FullNB->MakeImprint(theDetector->HallLog(), NoShiftR, &NoRotR, 0);

    }
    else if(shellStatus == "full" && !backwardShellStatus && !forwardShellStatus) {

       G4ThreeVector TarPos1
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos1.rotateY(ModuleEuler[0][1]);
      TarPos1.rotateZ(ModuleEuler[0][2]);
      G4RotationMatrix RotShell1 = G4RotationMatrix::IDENTITY;
      RotShell1.rotateZ(ModuleEuler[0][0]);
      RotShell1.rotateY( TarPos1.getTheta() );
      RotShell1.rotateZ( TarPos1.getPhi() );

      G4ThreeVector TarPos2
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos2.rotateY(ModuleEuler[1][1]);
      TarPos2.rotateZ(ModuleEuler[1][2]);
      G4RotationMatrix RotShell2 = G4RotationMatrix::IDENTITY;
      RotShell2.rotateZ(ModuleEuler[1][0]);
      RotShell2.rotateY( TarPos2.getTheta() );
      RotShell2.rotateZ( TarPos2.getPhi() );
      
      G4ThreeVector TarPos3
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos3.rotateY(ModuleEuler[2][1]);
      TarPos3.rotateZ(ModuleEuler[2][2]);
      G4RotationMatrix RotShell3 = G4RotationMatrix::IDENTITY;
      RotShell3.rotateZ(ModuleEuler[2][0]);
      RotShell3.rotateY( TarPos3.getTheta() );
      RotShell3.rotateZ( TarPos3.getPhi() );

      G4ThreeVector TarPos4
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos4.rotateY(ModuleEuler[3][1]);
      TarPos4.rotateZ(ModuleEuler[3][2]);
      G4RotationMatrix RotShell4 = G4RotationMatrix::IDENTITY;
      RotShell4.rotateZ(ModuleEuler[3][0]);
      RotShell4.rotateY( TarPos4.getTheta() );
      RotShell4.rotateZ( TarPos4.getPhi() );

      G4ThreeVector TarPos5
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos5.rotateY(ModuleEuler[4][1]);
      TarPos5.rotateZ(ModuleEuler[4][2]);
      G4RotationMatrix RotShell5 = G4RotationMatrix::IDENTITY;
      RotShell5.rotateZ(ModuleEuler[4][0]);
      RotShell5.rotateY( TarPos5.getTheta() );
      RotShell5.rotateZ( TarPos5.getPhi() );
      
      G4ThreeVector TarPos26
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos26.rotateY(ModuleEuler[25][1]);
      TarPos26.rotateZ(ModuleEuler[25][2]);
      G4RotationMatrix RotShell26 = G4RotationMatrix::IDENTITY;
      RotShell26.rotateZ(ModuleEuler[25][0]);
      RotShell26.rotateY( TarPos26.getTheta() );
      RotShell26.rotateZ( TarPos26.getPhi() );

      G4ThreeVector TarPos27
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos27.rotateY(ModuleEuler[26][1]);
      TarPos27.rotateZ(ModuleEuler[26][2]);
      G4RotationMatrix RotShell27 = G4RotationMatrix::IDENTITY;
      RotShell27.rotateZ(ModuleEuler[26][0]);
      RotShell27.rotateY( TarPos27.getTheta() );
      RotShell27.rotateZ( TarPos27.getPhi() );
      
      G4ThreeVector TarPos28
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos28.rotateY(ModuleEuler[27][1]);
      TarPos28.rotateZ(ModuleEuler[27][2]);
      G4RotationMatrix RotShell28 = G4RotationMatrix::IDENTITY;
      RotShell28.rotateZ(ModuleEuler[27][0]);
      RotShell28.rotateY( TarPos28.getTheta() );
      RotShell28.rotateZ( TarPos28.getPhi() );

      G4ThreeVector TarPos29
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos29.rotateY(ModuleEuler[28][1]);
      TarPos29.rotateZ(ModuleEuler[28][2]);
      G4RotationMatrix RotShell29 = G4RotationMatrix::IDENTITY;
      RotShell29.rotateZ(ModuleEuler[28][0]);
      RotShell29.rotateY( TarPos29.getTheta() );
      RotShell29.rotateZ( TarPos29.getPhi() );

      G4ThreeVector TarPos30
	= G4ThreeVector(0, 0, std::sqrt(std::pow(RBar,2)-std::pow(hPrime,2)));
      TarPos30.rotateY(ModuleEuler[29][1]);
      TarPos30.rotateZ(ModuleEuler[29][2]);
      G4RotationMatrix RotShell30 = G4RotationMatrix::IDENTITY;
      RotShell30.rotateZ(ModuleEuler[29][0]);
      RotShell30.rotateY( TarPos30.getTheta() );
      RotShell30.rotateZ( TarPos30.getPhi() );

      // Remove the forward central pentagon.
      G4SubtractionSolid* FCutPenta
	= new G4SubtractionSolid("FCutPenta", shellF, PentaCut,
				 G4Transform3D(NoRotR,
					       G4ThreeVector(0,0,halfheight)));
      // Remove the backward central pentagon.
      G4SubtractionSolid* FBCutPenta
	= new G4SubtractionSolid("FBCutPenta", FCutPenta, BPentaCut,
				 G4Transform3D(NoRotR,
					       G4ThreeVector(0,0,-halfheight)));
      
      // Make the forward flats.
      G4MultiUnion* FRingFlats = new G4MultiUnion("FRingFlats");
      G4Transform3D tr1 = G4Transform3D(RotShell1, TarPos1);
      FRingFlats->AddNode(*tripleShape, tr1);
      G4Transform3D tr2 = G4Transform3D(RotShell2, TarPos2);
      FRingFlats->AddNode(*tripleShape, tr2);
      G4Transform3D tr3 = G4Transform3D(RotShell3, TarPos3);
      FRingFlats->AddNode(*tripleShape, tr3);
      G4Transform3D tr4 = G4Transform3D(RotShell4, TarPos4);
      FRingFlats->AddNode(*tripleShape, tr4);
      G4Transform3D tr5 = G4Transform3D(RotShell5, TarPos5);
      FRingFlats->AddNode(*tripleShape, tr5);
      FRingFlats->Voxelize();

      // Cut the forward flats.
      G4SubtractionSolid* FBCutPentaFFlats
	= new G4SubtractionSolid("FRingCut", FBCutPenta, FRingFlats,
				 G4Transform3D(NoRotR,
					       G4ThreeVector(0,0,0)));
            
      // Make the backward flats.
      G4MultiUnion* BRingFlats = new G4MultiUnion("BRingFlats");
      
      G4Transform3D tr30 = G4Transform3D(RotShell30, TarPos30);
      BRingFlats->AddNode(*tripleShape, tr30);
      G4Transform3D tr29 = G4Transform3D(RotShell29, TarPos29);
      BRingFlats->AddNode(*tripleShape, tr29);
      G4Transform3D tr28 = G4Transform3D(RotShell28, TarPos28);
      BRingFlats->AddNode(*tripleShape, tr28);
      G4Transform3D tr27 = G4Transform3D(RotShell27, TarPos27);
      BRingFlats->AddNode(*tripleShape, tr27);
      G4Transform3D tr26 = G4Transform3D(RotShell26, TarPos26);
      BRingFlats->AddNode(*tripleShape, tr26);

      BRingFlats->Voxelize();

      // Cut the backward flats.
      G4SubtractionSolid* FBCutPentaFBFlats
      	= new G4SubtractionSolid("FBCutPentaFBFlats", FBCutPentaFFlats, BRingFlats,
      				 G4Transform3D(NoRotR,
      					       G4ThreeVector(0,0,0)));

      G4LogicalVolume* logicFBCutPentaFBFlats
	= new G4LogicalVolume(FBCutPentaFBFlats, matShell, "Shell_log", 0, 0, 0 );

      new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "MountingShell", logicFBCutPentaFBFlats,
      			theDetector->HallPhys(), false, 0 );

    }     
     
  }
  if(shellStatus == "full" && forwardShellStatus && backwardShellStatus) {
    //places the full Greta Shell with forward and backward rings
    G4LogicalVolume* logicShellF = new G4LogicalVolume(shellF, matShell,
						       "Shell_log", 0, 0, 0 );
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), "MountingShell", logicShellF,
		      theDetector->HallPhys(), false, 0 );
  }
  else if(shellStatus == "left" && forwardShellStatus && backwardShellStatus) {
    //places full Left Hemisphere with forward and backward rings
    G4AssemblyVolume* LeftHemi = new G4AssemblyVolume();
    LeftHemi->AddPlacedVolume(logicCutout30, NoShiftL, &NoRotL);
    LeftHemi->AddPlacedVolume(logicbump23, NoShiftL, &NoRotL);
    LeftHemi->AddPlacedVolume(logicbump3, NoShiftL, &NoRotL);
    LeftHemi->MakeImprint(theDetector->HallLog(), NoShiftL, &NoRotL, 0);
  }
  else if(shellStatus == "right" && forwardShellStatus && backwardShellStatus) {
    //places full Right Hemisphere with forward and backward rings
    G4AssemblyVolume* RightHemi = new G4AssemblyVolume();
    RightHemi->AddPlacedVolume(logicCutout3, NoShiftL, &NoRotL);
    RightHemi->AddPlacedVolume(logicbump10, NoShiftL, &NoRotL);
    RightHemi->AddPlacedVolume(logicbump30, NoShiftL, &NoRotL);
    RightHemi->MakeImprint(theDetector->HallLog(), NoShiftL, &NoRotL, 0);
  }
   
  if(shellStatus == "Greta")
    G4cout << "Constructed the GRETA mounting shell" << G4endl;
  else if(shellStatus == "Greta_Left")
    G4cout << "Constructed the left hemisphere of the GRETA mounting shell"
	   << G4endl;
  else if(shellStatus == "Greta_Right")
    G4cout << "Constructed the right hemisphere of the GRETA mounting shell"
	   << G4endl;
  if(!forwardShellStatus && !backwardShellStatus)
    G4cout << "   omitting the forward and backward rings"
	   << G4endl;
  else if(!forwardShellStatus)
    G4cout << "   omitting the forward ring"
	   << G4endl;
  else if(!backwardShellStatus)
    G4cout << "   omitting the backward ring"
	   << G4endl;
  
  G4cout << "  Shell radius: " << Rmin << " -- " << Rmax << G4endl;
  G4cout << "  Shell material: " << matShell->GetName() << G4endl;
  
}

//Creates and returns the full or partial Greta shell
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
    G4Exception("Greta_Shell::Shell()", "Error", FatalException,
		"half argument must be set to LEFT or RIGHT or FULL");

  G4ThreeVector Origin = G4ThreeVector(0, 0, 0);
  G4RotationMatrix NoRot = G4RotationMatrix::IDENTITY;
  
  G4Sphere* solidShell = new G4Sphere( "solidShell", Rmin, Rmax,
  				       Phi0, dPhi, 0., 180.*deg);

  G4Tubs* smallPort = new G4Tubs("smallPort", 0, smallPortRadius,
  				 1.5*Rmax,
  				 0., 360.*deg);

  // Beam ports
  G4SubtractionSolid* shell
    = new G4SubtractionSolid("Shell", solidShell, smallPort,
  			     G4Transform3D(Rot, Origin));

  // Other small ports
  
  // Opposites
  // 1 : 10
  // 2 :  9
  // 3 :  8
  // 4 :  7
  // 5 :  6

  G4int iMin = 0;
  G4int iMax = 5;

  for(G4int i = iMin; i<iMax; i++){
    if( abs(SmallPortStatus[i] <= 1 )) {
      Rot = G4RotationMatrix::IDENTITY;
      Rot.rotateY( PosSP[i].getTheta() );
      Rot.rotateZ( PosSP[i].getPhi() );
      shell = new G4SubtractionSolid ("Shell", shell, smallPort,
				      G4Transform3D(Rot, Origin));

    }
  }

  // Module ports
  G4Tubs* modulePort
    = new G4Tubs("modulePort", 0, modulePortRadius, 1.5*Rmax,
		   0., 360.*deg);

  // L 1,2,3S 6,7,     11,12,13,14,15, 21,22,23S, 26,27
  // R 4,5    8,9,10S, 16,17,18,19,20, 24,25,     28,29,30S

  // Opposites
  //  3 : 30
  //  2 : 29
  //  1 : 28
  //  4 : 26
  //  5 : 27
  
  //  6 : 24
  //  7 : 25
  //  8 : 21
  //  9 : 22
  // 10 : 23
  
  // 11 : 16
  // 12 : 17
  // 13 : 18
  // 14 : 19
  // 15 : 20

  iMin = 0;
  iMax = 15;

  for(G4int i = iMin; i<iMax; i++){
    if( ( ModulePortStatus[i] <= 0 ) || ( ModulePortStatus[i] >= 0 ) ) {
      Pos = Pos0;
      Pos.rotateZ(ModuleEuler[i][0]);
      Pos.rotateY(ModuleEuler[i][1]);
      Pos.rotateZ(ModuleEuler[i][2]);
      Rot = G4RotationMatrix::IDENTITY;
      Rot.rotateY( Pos.getTheta() );
      Rot.rotateZ( Pos.getPhi() );
      shell = new G4SubtractionSolid ("Shell", shell, modulePort,
				      G4Transform3D(Rot, Origin));
    }
  }
  
  return shell;

}
