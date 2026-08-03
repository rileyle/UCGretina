#include "Greta_Shell.hh"

Greta_Shell::Greta_Shell()
{

  // Dimensions taken from fabricationprint_25j2106a-1.pdf
  Rmin      = 1022.0/2.0*mm;
  Rmax      = 1276.0/2.0*mm;

  smallPortRadius  = 7.0/2.0*2.54*cm;    // DETAIL C, SHEET2
  modulePortRadius = 12.008/2.0*2.54*cm; // DETAIL A, SHEET2
  leftOffset       = 0.;
  rightOffset      = 0.;

  flatHeight = 540*mm;
  
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

  // Using the vertices surrounding Hole 2 from the CAD drawings and the
  // center of the forward port to define the flat shape. This can be used
  // to cut all 10 flats in the forward and backward rings.
  //
  // The same shape: (x, y, 540*mm) vectors, describe the
  // forward/backward ring cutouts/bumpouts (holes 3 and 30).
  //
  // (Geant4 x = Drawing y; Geant4 y = - Drawing x)
  CADPosFlats[0] = G4TwoVector(     0, -207.0); 
  CADPosFlats[1] = G4TwoVector( 203.2,  -79.8); 
  CADPosFlats[2] = G4TwoVector( 183.6,   92.7); 
  CADPosFlats[3] = G4TwoVector(     0,  210.0); 
  CADPosFlats[4] = G4TwoVector(-flatHeight*tan(31.717473*deg), 0); 

  // The corners of the polygon bounding holes 23 and 10 were used
  // to determine the shape of the irregular hexagon surrounding these
  // holes. This shape defines the extruded solid for generating cutouts
  // and bumpouts for holes 23 and 10 which are split by the hemisphere
  // boundary and have permanent cutouts/bumpouts on the hemispheres.
  // This is the upper surface of the extruded solid at
  // (Rmax + Rmin)/2 + Rmax - Rmin from the origin.
  // (Geant4 x = Drawing y; Geant4 y = - Drawing x)
  CADPosHex[0] = G4TwoVector(     0, -269.0); //      0, -266.0
  CADPosHex[1] = G4TwoVector(-235.8, -124.0); //  236.8, -119.0
  CADPosHex[2] = G4TwoVector(-260.6,  106.3); //  241.6,  114.3   (120.3)
  CADPosHex[3] = G4TwoVector(     0,  269.0); //      0,  266.0
  CADPosHex[4] = G4TwoVector( 235.8,  124.0); // -236.8,  119.0
  CADPosHex[5] = G4TwoVector( 260.6, -106.3); // -241.6, -114.3   (120.3)

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

//Modifies and places the GRETA mounting shell
void Greta_Shell::Placement(G4String shellStatus,
			    G4bool forwardShellStatus,
			    G4bool backwardShellStatus)
{

  G4ThreeVector ShiftL = G4ThreeVector(leftOffset, 0, 0);
  G4ThreeVector ShiftR = G4ThreeVector(-rightOffset, 0, 0);
  G4RotationMatrix NoRot = G4RotationMatrix::IDENTITY;

  G4RunManager* runManager = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector
    = (DetectorConstruction*) runManager->GetUserDetectorConstruction();

  if(shellStatus == "full") {
    Hemi("left",
	 forwardShellStatus,
	 backwardShellStatus)->MakeImprint(theDetector->HallLog(),
					   ShiftL, &NoRot, 0);
    Hemi("right",
	 forwardShellStatus,
	 backwardShellStatus)->MakeImprint(theDetector->HallLog(),
					   ShiftR, &NoRot, 0);
    G4cout << "Constructed the GRETA mounting shell" << G4endl;
    G4cout << "   left hemisphere offset " << leftOffset*mm << " mm"
	   << G4endl;
    G4cout << "   right hemisphere offset " << rightOffset*mm << " mm"
	   << G4endl;
  } else if(shellStatus == "left"){
    Hemi("left",
	 forwardShellStatus,
	 backwardShellStatus)->MakeImprint(theDetector->HallLog(),
					   ShiftL, &NoRot, 0);
    G4cout << "Constructed the left hemisphere of the GRETA mounting shell"
	   << G4endl;
    G4cout << "   left hemisphere offset " << leftOffset*mm << " mm"
	   << G4endl;
  } else if(shellStatus == "right"){
    Hemi("right",
	 forwardShellStatus,
	 backwardShellStatus)->MakeImprint(theDetector->HallLog(),
					   ShiftR, &NoRot, 0);
    G4cout << "Constructed the right hemisphere of the GRETA mounting shell"
	   << G4endl;
    G4cout << "   right hemisphere offset " << rightOffset*mm << " mm"
	   << G4endl;
  }

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

// Creates and returns a sphere or hemisphere with module ports
// and small ports.
G4SubtractionSolid* Greta_Shell::Shell(G4String half)
{
  G4double Phi0=0;
  G4double dPhi=0;
  Rot = G4RotationMatrix::IDENTITY;

  if( FindMaterials() ) return NULL;
 
  if(half == "left"){
    Phi0 = -90.*deg;
    dPhi = 180.*deg;
  }
  else if(half == "right"){
    Phi0 =  90.*deg;
    dPhi = 180.*deg;
  }
  else if(half == "full"){
    Phi0 = 0.*deg;
    dPhi = 360.*deg;
  } else
    G4Exception("Greta_Shell::Shell()", "Error", FatalException,
		"half argument must be set to left or right or full");

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

// Builds the left or right hemisphere with cutouts and bumpouts
// and optional removal of the forward and/or backward ring.
G4AssemblyVolume* Greta_Shell::Hemi(G4String half,
				    G4bool forwardShellStatus,
				    G4bool backwardShellStatus)
{
  // Half the height of the extruded solids that cut fully through the shell.
  // We want these to extend from within Rmin to beyond Rmax. The solids have
  // twice the thickness of the shell.
  G4double halfheight = Rmax - Rmin;

  G4ThreeVector NoShift = G4ThreeVector(0, 0, 0);
  G4RotationMatrix NoRot = G4RotationMatrix::IDENTITY;
  
  std::vector<G4int> fFlats; // stores slot numbers of forward flats to cut
  std::vector<G4int> bFlats; // stores slot numbers of backward flats to cut
  G4int bumps[2] = {-1, -1};
  G4int cuts[2]  = {-1, -1};
  
  if(half == "left"){
    bumps[0] = 22;
    bumps[1] = 2;
    cuts[0]  = 9;
    cuts[1]  = 29;
    if(!forwardShellStatus){
      fFlats.push_back(0);
      fFlats.push_back(1);
      fFlats.push_back(2);
    }
    if(!backwardShellStatus){
      bFlats.push_back(25);
      bFlats.push_back(26);
      bFlats.push_back(29);
    }
  } else if(half == "right"){
    bumps[0] = 9;
    bumps[1] = 29;
    cuts[0]  = 22;
    cuts[1]  = 2;
  if(!forwardShellStatus){
      fFlats.push_back(2);
      fFlats.push_back(3);
      fFlats.push_back(4);
    }
    if(!backwardShellStatus){
      bFlats.push_back(27);
      bFlats.push_back(28);
      bFlats.push_back(29);
    }
  }
  else
    G4cout << "Greta_Shell::Hemi half parameter " << half << " is not defined "
	   << "(expecting left or right)."
	   << G4endl;

  G4AssemblyVolume* Hemisphere = new G4AssemblyVolume();
  
  G4SubtractionSolid *shell = Shell(half);

  // We define flatPoly in this scope, because it describes both flats and
  // forward/backward cutouts and bumpouts (holes 3 and 30).
  std::vector<G4TwoVector> flatPoly(5);
  flatPoly[0] = G4TwoVector(CADPosFlats[0].x(), CADPosFlats[0].y());
  flatPoly[1] = G4TwoVector(CADPosFlats[1].x(), CADPosFlats[1].y());
  flatPoly[2] = G4TwoVector(CADPosFlats[2].x(), CADPosFlats[2].y());
  flatPoly[3] = G4TwoVector(CADPosFlats[3].x(), CADPosFlats[3].y());
  flatPoly[4] = G4TwoVector(CADPosFlats[4].x(), CADPosFlats[4].y());
  
  // Remove forward / backward rings as requested.
  if(!forwardShellStatus || !backwardShellStatus){
    G4double innercut = 751*mm;
    if(!forwardShellStatus){
      // Subtract forward pentagon
      std::vector<G4TwoVector> pentagon(5);
      for(G4int i=0; i<5; i++)
	pentagon[i]
	  = G4TwoVector(innercut*sin(ModuleEuler[i][1])*cos(ModuleEuler[i][2]),
			innercut*sin(ModuleEuler[i][1])*sin(ModuleEuler[i][2]));
      G4ExtrudedSolid* PentaCut
	= new G4ExtrudedSolid("PentaCut", pentagon, halfheight,
			      G4TwoVector(0, 0),
			      ((Rmin+Rmax)/2 - halfheight)/((Rmin+Rmax)/2 + halfheight),
			      G4TwoVector(0, 0), 1);
      shell = new G4SubtractionSolid("shell", shell, PentaCut,
				     G4Transform3D(NoRot,
						   G4ThreeVector(0, 0,
								 (Rmin+Rmax)/2)));
    }
    if(!backwardShellStatus){
      // Subtract backward pentagon
      std::vector<G4TwoVector> pentagon(5);
      for(G4int i=0; i<5; i++)
	pentagon[i]
	  = G4TwoVector(innercut*sin(ModuleEuler[i+25][1])*cos(ModuleEuler[i+25][2]),
			innercut*sin(ModuleEuler[i+25][1])*sin(ModuleEuler[i+25][2]));
      G4RotationMatrix RotBackPenta = G4RotationMatrix::IDENTITY;
      RotBackPenta.rotateY(180*degree);
      G4ExtrudedSolid* PentaCut
	= new G4ExtrudedSolid("PentaCut", pentagon, halfheight,
			      G4TwoVector(0, 0),
			      ((Rmin+Rmax)/2 - halfheight)/((Rmin+Rmax)/2 + halfheight),
			      G4TwoVector(0, 0), 1);
      shell = new G4SubtractionSolid("shell", shell, PentaCut,
				     G4Transform3D(RotBackPenta,
						   G4ThreeVector(0, 0,
								 -(Rmin+Rmax)/2)));
    }

    // Subtract the flats
    std::vector<G4ExtrudedSolid::ZSection> zsections;
    zsections.push_back(G4ExtrudedSolid::ZSection(0,
						  G4TwoVector(0,0), 1));
    zsections.push_back(G4ExtrudedSolid::ZSection(Rmax - Rmin,
						  G4TwoVector(0,0),
						  (flatHeight + halfheight)/flatHeight));
    G4ExtrudedSolid* flatShape
      = new G4ExtrudedSolid("flatShape", flatPoly, zsections);

    std::vector<G4ExtrudedSolid::ZSection> zsections2;
    zsections2.push_back(G4ExtrudedSolid::ZSection(-halfheight,
						   G4TwoVector(0,0),
						   (flatHeight - halfheight)/flatHeight));
    zsections2.push_back(G4ExtrudedSolid::ZSection(halfheight,
						   G4TwoVector(0,0),
						   (flatHeight + halfheight)/flatHeight));
    G4ExtrudedSolid* doubleHoleShape
      = new G4ExtrudedSolid("doubleHoleShape", flatPoly, zsections2);
    
    if(!forwardShellStatus){
      G4MultiUnion* forwardFlats = new G4MultiUnion("ForwardFlats");
      for(auto i: fFlats){
	G4ThreeVector PosFlat = G4ThreeVector(0, 0, flatHeight);
	PosFlat.rotateY(ModuleEuler[i][1]);
	PosFlat.rotateZ(ModuleEuler[i][2]);
	G4RotationMatrix RotFlat = G4RotationMatrix::IDENTITY;
	RotFlat.rotateZ(ModuleEuler[i][0]);
	RotFlat.rotateY( PosFlat.getTheta() );
	RotFlat.rotateZ( PosFlat.getPhi() );
	G4Transform3D tr = G4Transform3D(RotFlat, PosFlat);
	if(half == "right" && i == 2)
	  forwardFlats->AddNode(*doubleHoleShape, tr);
	else
	  forwardFlats->AddNode(*flatShape, tr);
      }
      forwardFlats->Voxelize();
      shell
	= new G4SubtractionSolid("shell", shell, forwardFlats,
				 G4Transform3D(NoRot, NoShift));
    }
    if(!backwardShellStatus){
      G4MultiUnion* backwardFlats = new G4MultiUnion("BackwardFlats");
      for(auto i: bFlats){
	G4ThreeVector PosFlat = G4ThreeVector(0, 0, flatHeight);
	PosFlat.rotateY(ModuleEuler[i][1]);
	PosFlat.rotateZ(ModuleEuler[i][2]);
	G4RotationMatrix RotFlat = G4RotationMatrix::IDENTITY;
	RotFlat.rotateZ(ModuleEuler[i][0]);
	RotFlat.rotateY( PosFlat.getTheta() );
	RotFlat.rotateZ( PosFlat.getPhi() );
	G4Transform3D tr = G4Transform3D(RotFlat, PosFlat);
	if(half == "left" && i == 29)
	  backwardFlats->AddNode(*doubleHoleShape, tr);
	else
	  backwardFlats->AddNode(*flatShape, tr);
      }
      backwardFlats->Voxelize();
      shell
	= new G4SubtractionSolid("shell", shell, backwardFlats,
				 G4Transform3D(NoRot, NoShift));
    }
    
  }

  // Bumpouts and cutouts
  G4IntersectionSolid* bump1, *bump2;
  G4LogicalVolume* logicBump1, *logicBump2; 
  if(half == "left" || half == "right"){
    std::vector<G4TwoVector> hexagon(6);
    hexagon[0] = G4TwoVector(CADPosHex[0].x(), CADPosHex[0].y());
    hexagon[1] = G4TwoVector(CADPosHex[1].x(), CADPosHex[1].y());
    hexagon[2] = G4TwoVector(CADPosHex[2].x(), CADPosHex[2].y());
    hexagon[3] = G4TwoVector(CADPosHex[3].x(), CADPosHex[3].y());
    hexagon[4] = G4TwoVector(CADPosHex[4].x(), CADPosHex[4].y());
    hexagon[5] = G4TwoVector(CADPosHex[5].x(), CADPosHex[5].y());

    // Middle cutout
    G4ExtrudedSolid* extrudedHexagon
      = new G4ExtrudedSolid("extrudedHexagon", hexagon, halfheight,
			    G4TwoVector(0, 0),
			    ((Rmin+Rmax)/2 - halfheight)/((Rmin+Rmax)/2 + halfheight),
			    G4TwoVector(0, 0), 1);

    G4ThreeVector PosHex = G4ThreeVector(0, 0, (Rmin+Rmax)/2);
    PosHex.rotateZ(ModuleEuler[cuts[0]][0]);
    PosHex.rotateY(ModuleEuler[cuts[0]][1]);
    PosHex.rotateZ(ModuleEuler[cuts[0]][2]);
    G4RotationMatrix RotHex = G4RotationMatrix::IDENTITY;
    RotHex.rotateZ( ModuleEuler[cuts[0]][0] );
    RotHex.rotateY( PosHex.getTheta() );
    RotHex.rotateZ( PosHex.getPhi() );
    shell = new G4SubtractionSolid("shell", shell, extrudedHexagon,
				   G4Transform3D(RotHex, PosHex));

    // Middle bumpout
    G4SubtractionSolid* otherShell;
    if(half == "left")
      otherShell = Shell("right");
    else
      otherShell = Shell("left");

    G4ThreeVector PosBump1 = G4ThreeVector(0, 0, (Rmin+Rmax)/2);
    PosBump1.rotateZ(ModuleEuler[bumps[0]][0]);
    PosBump1.rotateY(ModuleEuler[bumps[0]][1]);
    PosBump1.rotateZ(ModuleEuler[bumps[0]][2]);
    G4RotationMatrix RotBump1 = G4RotationMatrix::IDENTITY;
    RotBump1.rotateZ( ModuleEuler[bumps[0]][0] );
    RotBump1.rotateY( PosBump1.getTheta() );
    RotBump1.rotateZ( PosBump1.getPhi() );

    bump1 = new G4IntersectionSolid("bump1", otherShell, extrudedHexagon,
				    G4Transform3D(RotBump1, PosBump1));
    logicBump1 = new G4LogicalVolume(bump1, matShell, "Bump1_log", 0, 0, 0 );

    Hemisphere->AddPlacedVolume(logicBump1, NoShift, &NoRot);
    
    // Forward/backward cutout
    G4ExtrudedSolid* extrudedPentagon
      = new G4ExtrudedSolid("extrudedPentagon", flatPoly, flatHeight,
			    G4TwoVector(0, 0), (flatHeight - halfheight)/flatHeight,
			    G4TwoVector(0, 0), (flatHeight + halfheight)/flatHeight);

    G4ThreeVector PosPent = G4ThreeVector(0, 0, flatHeight);
    PosPent.rotateZ(ModuleEuler[cuts[1]][0]);
    PosPent.rotateY(ModuleEuler[cuts[1]][1]);
    PosPent.rotateZ(ModuleEuler[cuts[1]][2]);
    G4RotationMatrix RotPent = G4RotationMatrix::IDENTITY;
    RotPent.rotateZ( ModuleEuler[cuts[1]][0] );
    RotPent.rotateY( PosPent.getTheta() );
    RotPent.rotateZ( PosPent.getPhi() );
    shell = new G4SubtractionSolid("shell", shell, extrudedPentagon,
				   G4Transform3D(RotPent, PosPent));

    // Forward/backward bumpout
    if( (bumps[1] <= 4  &&  forwardShellStatus) ||
	(bumps[1] >= 25 && backwardShellStatus) ){
      G4ThreeVector PosBump2 = G4ThreeVector(0, 0, flatHeight);
      PosBump2.rotateZ(ModuleEuler[bumps[1]][0]);
      PosBump2.rotateY(ModuleEuler[bumps[1]][1]);
      PosBump2.rotateZ(ModuleEuler[bumps[1]][2]);
      G4RotationMatrix RotBump2 = G4RotationMatrix::IDENTITY;
      RotBump2.rotateZ( ModuleEuler[bumps[1]][0] );
      RotBump2.rotateY( PosBump2.getTheta() );
      RotBump2.rotateZ( PosBump2.getPhi() );

      bump2 = new G4IntersectionSolid("bump2", otherShell, extrudedPentagon,
				      G4Transform3D(RotBump2, PosBump2));
      logicBump2 = new G4LogicalVolume(bump2, matShell, "Bump1_log", 0, 0, 0 );

      Hemisphere->AddPlacedVolume(logicBump2, NoShift, &NoRot);
    }
  }
  
  G4LogicalVolume* logicShell
    = new G4LogicalVolume(shell, matShell, "Shell_log", 0, 0, 0 );

  Hemisphere->AddPlacedVolume(logicShell, NoShift, &NoRot);

  return Hemisphere;

}
