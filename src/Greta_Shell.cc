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

  // Module Port Euler angles (relative to Slot 0)
  // Psi                              Slot   Hemisphere
  // Theta
  // Phi
  
  ModuleEuler[0][0] =  288.00*deg;    //  0    North
  ModuleEuler[0][1] =    0.00*deg;
  ModuleEuler[0][2] =    0.00*deg;
  ModuleEuler[1][0] =    0.00*deg;    //  1    North
  ModuleEuler[1][1] =    0.00*deg;
  ModuleEuler[1][2] =    0.00*deg;
  ModuleEuler[2][0] =   72.00*deg;    //  2    Split
  ModuleEuler[2][1] =    0.00*deg;
  ModuleEuler[2][2] =    0.00*deg;
  ModuleEuler[3][0] =  144.00*deg;    //  3    South
  ModuleEuler[3][1] =    0.00*deg;
  ModuleEuler[3][2] =    0.00*deg;
  ModuleEuler[4][0] =  216.00*deg;    //  4    South
  ModuleEuler[4][1] =    0.00*deg;
  ModuleEuler[4][2] =    0.00*deg;
  ModuleEuler[5][0] =  235.58*deg;    //  5    North
  ModuleEuler[5][1] =   63.43495*deg;
  ModuleEuler[5][2] =   16.420*deg;
  ModuleEuler[6][0] =   91.58*deg;    //  6    North
  ModuleEuler[6][1] =   63.43495*deg;
  ModuleEuler[6][2] =   16.420*deg;
  ModuleEuler[7][0] =   91.58*deg;    //  7    South
  ModuleEuler[7][1] =   63.43495*deg;
  ModuleEuler[7][2] =   88.420*deg;
  ModuleEuler[8][0] =   91.58*deg;    //  8    South
  ModuleEuler[8][1] =   63.43495*deg;
  ModuleEuler[8][2] =  160.420*deg;
  ModuleEuler[9][0] =   91.58*deg;    //  9    Split
  ModuleEuler[9][1] =   63.43495*deg;
  ModuleEuler[9][2] =  232.420*deg;
  ModuleEuler[10][0] = -52.42*deg;    // 10    North
  ModuleEuler[10][1] =   63.43495*deg;
  ModuleEuler[10][2] =  304.420*deg;
  ModuleEuler[11][0] = -340.42*deg;   // 11    North
  ModuleEuler[11][1] =   63.43495*deg;
  ModuleEuler[11][2] =  304.420*deg;
  ModuleEuler[12][0] =  -52.42*deg;   // 12    North
  ModuleEuler[12][1] =   63.43495*deg;
  ModuleEuler[12][2] =   16.420*deg;
  ModuleEuler[13][0] =   19.58*deg;   // 13    North
  ModuleEuler[13][1] =   63.43495*deg;
  ModuleEuler[13][2] =   16.420*deg;
  ModuleEuler[14][0] =  -52.42*deg;   // 14    North
  ModuleEuler[14][1] =   63.43495*deg;
  ModuleEuler[14][2] =   88.420*deg;
  ModuleEuler[15][0] =   19.58*deg;   // 15    South
  ModuleEuler[15][1] =   63.43495*deg;
  ModuleEuler[15][2] =   88.420*deg;
  ModuleEuler[16][0] =   -52.42*deg;  // 16    South
  ModuleEuler[16][1] =   63.43495*deg;
  ModuleEuler[16][2] =   160.420*deg;
  ModuleEuler[17][0] =   19.58*deg;   // 17    South
  ModuleEuler[17][1] =   63.43495*deg;
  ModuleEuler[17][2] =   160.420*deg;
  ModuleEuler[18][0] =   -52.42*deg;  // 18    South
  ModuleEuler[18][1] =   63.43495*deg;
  ModuleEuler[18][2] =   232.420*deg;
  ModuleEuler[19][0] =   19.58*deg;   // 19    South
  ModuleEuler[19][1] =   63.43495*deg;
  ModuleEuler[19][2] =   232.420*deg;
  ModuleEuler[20][0] =   55.58*deg;   // 20    North
  ModuleEuler[20][1] =  116.56505*deg;
  ModuleEuler[20][2] =  268.420*deg;
  ModuleEuler[21][0] =  -88.42*deg;   // 21    North
  ModuleEuler[21][1] =  116.56505*deg;
  ModuleEuler[21][2] =   52.420*deg;
  ModuleEuler[22][0] =   55.58*deg;   // 22    Split
  ModuleEuler[22][1] =  116.56505*deg;
  ModuleEuler[22][2] =   52.420*deg;
  ModuleEuler[23][0] =   55.58*deg;   // 23    South
  ModuleEuler[23][1] =  116.56505*deg;
  ModuleEuler[23][2] =  124.420*deg;
  ModuleEuler[24][0] =   55.58*deg;   // 24    South
  ModuleEuler[24][1] =  116.56505*deg;
  ModuleEuler[24][2] =  196.420*deg;
  ModuleEuler[25][0] =  -16.42*deg;   // 25    North
  ModuleEuler[25][1] =  116.56505*deg;
  ModuleEuler[25][2] =  340.420*deg;
  ModuleEuler[26][0] =  -16.42*deg;   // 26    North
  ModuleEuler[26][1] =  116.56505*deg;
  ModuleEuler[26][2] =   52.420*deg;
  ModuleEuler[27][0] =   -16.42*deg;  // 27    South
  ModuleEuler[27][1] =   116.56505*deg;
  ModuleEuler[27][2] =   124.420*deg;
  ModuleEuler[28][0] =   -16.42*deg;  // 28    South
  ModuleEuler[28][1] =   116.56505*deg;
  ModuleEuler[28][2] =   196.420*deg;
  ModuleEuler[29][0] =  -16.42*deg;   // 29    Split
  ModuleEuler[29][1] =  116.56505*deg;
  ModuleEuler[29][2] =  268.420*deg;

  matShell       = NULL;
  matShellName   = "Al";

  // This is the postion of the center of the front faces of a cluster before
  // it is placed according to its euler angles (the center of Slot 0).
  Pos0.setX(90.734*mm);
  Pos0.setY(26.7389*mm);
  Pos0.setZ(153.053*mm);
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
       status != "GretaLH" ){
    G4cout << "Shell status " << status << " is not defined." << G4endl;
    return;
  }

  if( FindMaterials() ) return;

  G4RunManager* runManager               = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector = (DetectorConstruction*) runManager->GetUserDetectorConstruction();
  
  //////////////////////////////////////////////////
  // The GRETA mounting shell
  //////////////////////////////////////////////////

  G4double Phi0 =   0.*deg;
  G4double dPhi = 360.*deg;

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
    shell = new G4SubtractionSolid ("Shell",shell,smallPort,G4Transform3D(Rot,PosSP[i]));
  }

  // Make space for the LH target
  if( status == "GretaLH"){
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

    shell = new G4SubtractionSolid ("Shell",shell,modulePort,G4Transform3D(Rot,Pos));
  }

  Rot = G4RotationMatrix::IDENTITY;

  G4LogicalVolume* logicShell = new G4LogicalVolume(shell, matShell, "Shell_log", 0, 0, 0 );

  new G4PVPlacement(0, G4ThreeVector(), "mountingShell", logicShell, theDetector->HallPhys(), false, 0 );

  G4VisAttributes *pVA = new G4VisAttributes( G4Colour(0.0, 1.0, 1.0) );
  logicShell->SetVisAttributes( pVA );

  G4cout << "Constructed the " << status << " shell." << G4endl;
  G4cout << "  Shell radius: " << Rmin << " -- " << Rmax << G4endl;
  G4cout << "  Shell material: " << matShell->GetName() << G4endl;

}
