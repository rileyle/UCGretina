#include "Greta_Shell.hh"

Greta_Shell::Greta_Shell()
{

  // Dimensions taken from fabricationprint_25j2106a-1.pdf
  Rmin      = 1022.0/2.0*mm;
  Rmax      = 1276.0/2.0*mm;

  //  mountPortRadius   = 3.0/2.0*2.54*cm;   // DETAIL F, SHEET2
  smallPortRadius   = 7.0/2.0*2.54*cm;   // DETAIL C, SHEET2
  modulePortRadius = 12.008/2.0*2.54*cm; // DETAIL A, SHEET2

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

  // North: -1,  Split: 0,  South: 1
  SmallPortStatus[0] = -1;
  SmallPortStatus[1] = 1;
  SmallPortStatus[2] = 1;
  SmallPortStatus[3] = -1;
  SmallPortStatus[4] = -1;
  SmallPortStatus[5] = 1;
  SmallPortStatus[6] = 1;
  SmallPortStatus[7] = -1;
  SmallPortStatus[8] = -1;
  SmallPortStatus[9] = 1;

  // Module Port Euler angles (relative to Slot 0)
  // Psi                              Slot   Hemisphere
  // Theta
  // Phi
  
  ModuleEuler[0][0] =     0.000000*deg;   //  0    North
  ModuleEuler[0][1] =    31.717473*deg;
  ModuleEuler[0][2] =   -54.000000*deg;
  ModuleEuler[1][0] =     0.000000*deg;   //  1    North
  ModuleEuler[1][1] =    31.717473*deg;
  ModuleEuler[1][2] =    18.000000*deg;
  ModuleEuler[2][0] =     0.000000*deg;   //  2    Split
  ModuleEuler[2][1] =    31.717473*deg;
  ModuleEuler[2][2] =    90.000000*deg;
  ModuleEuler[3][0] =     0.000000*deg;   //  3    South
  ModuleEuler[3][1] =    31.717473*deg;
  ModuleEuler[3][2] =   162.000000*deg;
  ModuleEuler[4][0] =     0.000000*deg;   //  4    South
  ModuleEuler[4][1] =    31.717473*deg;
  ModuleEuler[4][2] =  -126.000000*deg;
  ModuleEuler[5][0] =    90.000000*deg;   //  5    North
  ModuleEuler[5][1] =    58.282526*deg;
  ModuleEuler[5][2] =   -18.000000*deg;
  ModuleEuler[6][0] =    90.000000*deg;   //  6    North
  ModuleEuler[6][1] =    58.282526*deg;
  ModuleEuler[6][2] =    54.000000*deg;
  ModuleEuler[7][0] =    90.000000*deg;   //  7    South
  ModuleEuler[7][1] =    58.282526*deg;
  ModuleEuler[7][2] =   126.000000*deg;
  ModuleEuler[8][0] =    90.000000*deg;   //  8    South
  ModuleEuler[8][1] =    58.282526*deg;
  ModuleEuler[8][2] =  -162.000000*deg;
  ModuleEuler[9][0] =    90.000000*deg;   //  9    Split
  ModuleEuler[9][1] =    58.282526*deg;
  ModuleEuler[9][2] =   -90.000000*deg;
  ModuleEuler[10][0] =  -31.750000*deg;   // 10    North
  ModuleEuler[10][1] =   90.000000*deg;
  ModuleEuler[10][2] =  -72.000000*deg;
  ModuleEuler[11][0] =   31.750000*deg;   // 11    North
  ModuleEuler[11][1] =   90.000000*deg;
  ModuleEuler[11][2] =  -36.000000*deg;
  ModuleEuler[12][0] =  -31.750000*deg;   // 12    North
  ModuleEuler[12][1] =   90.000000*deg;
  ModuleEuler[12][2] =    0.000000*deg;
  ModuleEuler[13][0] =   31.750000*deg;   // 13    North
  ModuleEuler[13][1] =   90.000000*deg;
  ModuleEuler[13][2] =   36.000000*deg;
  ModuleEuler[14][0] =  -31.750000*deg;   // 14    North
  ModuleEuler[14][1] =   90.000000*deg;
  ModuleEuler[14][2] =   72.000000*deg;
  ModuleEuler[15][0] =   31.750000*deg;   // 15    South
  ModuleEuler[15][1] =   90.000000*deg;
  ModuleEuler[15][2] =  108.000000*deg;
  ModuleEuler[16][0] =  -31.750000*deg;   // 16    South
  ModuleEuler[16][1] =   90.000000*deg;
  ModuleEuler[16][2] =  144.000000*deg;
  ModuleEuler[17][0] =   31.750000*deg;   // 17    South
  ModuleEuler[17][1] =   90.000000*deg;
  ModuleEuler[17][2] =  180.000000*deg;
  ModuleEuler[18][0] =  -31.750000*deg;   // 18    South
  ModuleEuler[18][1] =   90.000000*deg;
  ModuleEuler[18][2] = -144.000000*deg;
  ModuleEuler[19][0] =   31.750000*deg;   // 19    South
  ModuleEuler[19][1] =   90.000000*deg; 
  ModuleEuler[19][2] = -108.000000*deg;
  ModuleEuler[20][0] =  -90.000000*deg;   // 20    North
  ModuleEuler[20][1] =  121.717475*deg;
  ModuleEuler[20][2] =  -54.000000*deg;
  ModuleEuler[21][0] =  -90.000000*deg;   // 21    North
  ModuleEuler[21][1] =  121.717475*deg;
  ModuleEuler[21][2] =   18.000000*deg;
  ModuleEuler[22][0] =  -90.000000*deg;   // 22    Split
  ModuleEuler[22][1] =  121.717475*deg;
  ModuleEuler[22][2] =   90.000000*deg;
  ModuleEuler[23][0] =  -90.000000*deg;   // 23    South
  ModuleEuler[23][1] =  121.717475*deg;
  ModuleEuler[23][2] =  162.000000*deg;
  ModuleEuler[24][0] =  -90.000000*deg;   // 24    South
  ModuleEuler[24][1] =  121.717475*deg;
  ModuleEuler[24][2] = -126.000000*deg;
  ModuleEuler[25][0] =  180.000000*deg;   // 25    North
  ModuleEuler[25][1] =  148.282526*deg;
  ModuleEuler[25][2] =  -18.000000*deg;
  ModuleEuler[26][0] =  180.000000*deg;   // 26    North
  ModuleEuler[26][1] =  148.282526*deg;
  ModuleEuler[26][2] =   54.000000*deg;
  ModuleEuler[27][0] =  180.000000*deg;   // 27    South
  ModuleEuler[27][1] =  148.282526*deg;
  ModuleEuler[27][2] =  126.000000*deg;
  ModuleEuler[28][0] =  180.000000*deg;   // 28    South
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
  //  Pos0.setX(90.734*mm);
  //  Pos0.setY(26.7389*mm);
  //  Pos0.setZ(153.053*mm);
  // Rescale to place the Tubs at a radial distance from the origin in the 
  // middle of the mounting shell. 
  Pos0.setMag((Rmin+Rmax)/2.*mm);

  Rot0 = G4RotationMatrix::IDENTITY;
  // This orients the cluster before it is placed according to its 
  // euler angles.
  Rot0.rotateY( Pos0.getTheta() );
  Rot0.rotateZ( Pos0.getPhi() );

  Rot = G4RotationMatrix::IDENTITY;

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

  if( FindMaterials() ) return;

  G4RunManager* runManager               = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector = (DetectorConstruction*) runManager->GetUserDetectorConstruction();
  
  //////////////////////////////////////////////////
  // The GRETA mounting shell
  //////////////////////////////////////////////////

  G4double Phi0 = 0*deg;
  G4double dPhi = 360*deg;

  if (status == "Greta_North" || status == "GretaLH_North" ) {
	Phi0 =  -90.*deg;
        dPhi = 180.*deg;
  }

  if (status == "Greta_South" || status == "GretaLH_South" ) {
	Phi0 =  90.*deg;
        dPhi = 180.*deg;
  }

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
    if( ( ( status == "Greta_North" || status == "GretaLH_North" ) && ModulePortStatus[i] == -1 ) ||
        ( ( status == "Greta_South" || status == "GretaLH_South" ) && ModulePortStatus[i] ==  1 ) ||
	status == "Greta" ||
        status == "GretaLH" )
      shell = new G4SubtractionSolid ("Shell",shell,smallPort,G4Transform3D(Rot,PosSP[i]));
  } 

  // Make space for the LH target
  if( status == "GretaLH" || status == "GretaLH_North" || status == "GretaLH_South"){
    G4Box* PortBox = new G4Box("PortBox", modulePortRadius, Rmax-Rmin, modulePortRadius);

    Rot = G4RotationMatrix::IDENTITY;

    shell = new G4SubtractionSolid ("Shell",shell,PortBox,
				    G4Transform3D(Rot,
						  G4ThreeVector(0.,
								(Rmax+Rmin)/2.,
								0.0)));
  }

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

    if( ( ( status == "Greta_North" || status == "GretaLH_North" )
	  && ModulePortStatus[i] == -1 ) ||
        ( ( status == "Greta_South" || status == "GretaLH_South" ) 
	  && ModulePortStatus[i] ==  1 ) ||
	status == "Greta" ||
        status == "GretaLH" )
          shell = new G4SubtractionSolid ("Shell",shell,modulePort,G4Transform3D(Rot,Pos));	

  }

  G4LogicalVolume* logicShell = new G4LogicalVolume(shell, matShell, "Shell_log", 0, 0, 0 );

  new G4PVPlacement(0, G4ThreeVector(), "mountingShell", logicShell, theDetector->HallPhys(), false, 0 );

  G4VisAttributes *pVA = new G4VisAttributes( G4Colour(0.0, 1.0, 1.0) );
  logicShell->SetVisAttributes( pVA );

  G4cout << "Constructed the " << status << " shell." << G4endl;
  G4cout << "  Shell radius: " << Rmin << " -- " << Rmax << G4endl;
  G4cout << "  Shell material: " << matShell->GetName() << G4endl;

}
