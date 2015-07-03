#include "Gretina_NSCL_Shell.hh"

Gretina_NSCL_Shell::Gretina_NSCL_Shell()
{

  // Dimensions taken from fabricationprint_25j2106a-1.pdf
  Rmin      = 1022.0/2.0*mm;
  Rmax      = 1276.0/2.0*mm;

  mountPortRadius   = 3.0/2.0*2.54*cm;   // DETAIL F, SHEET2
  smallPortRadius   = 7.0/2.0*2.54*cm;   // DETAIL C, SHEET2
  modulePortRadius = 12.008/2.0*2.54*cm; // DETAIL A, SHEET2
  S800PortRadius   = 13.287*2.54*cm;     // SHEET 1
  notchWidth       = 15.59*2.54*cm;
  notchThickness   = 4.921*2.54*cm+10*cm;

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
  // Psi                                   Slot   Hemisphere
  // Theta
  // Phi
  ModuleEuler[0][0] =    90.000000*deg;   //  7    South
  ModuleEuler[0][1] =    58.282526*deg;
  ModuleEuler[0][2] =   126.000000*deg;
  ModuleEuler[1][0] =    90.000000*deg;   //  8    South
  ModuleEuler[1][1] =    58.282526*deg;
  ModuleEuler[1][2] =  -162.000000*deg;
  ModuleEuler[2][0] =    31.750000*deg;   // 15    South
  ModuleEuler[2][1] =    90.000000*deg;
  ModuleEuler[2][2] =   108.000000*deg;
  ModuleEuler[3][0] =   -31.750000*deg;   // 16    South
  ModuleEuler[3][1] =    90.000000*deg;
  ModuleEuler[3][2] =   144.000000*deg;
  ModuleEuler[4][0] =   -31.750000*deg;   // 18    South
  ModuleEuler[4][1] =    90.000000*deg;
  ModuleEuler[4][2] =  -144.000000*deg;
  ModuleEuler[5][0] =    31.750000*deg;   // 19    South
  ModuleEuler[5][1] =    90.000000*deg; 
  ModuleEuler[5][2] =  -108.000000*deg;
  ModuleEuler[6][0] =   -90.000000*deg;   // 23    South
  ModuleEuler[6][1] =   121.717475*deg;
  ModuleEuler[6][2] =   162.000000*deg;
  ModuleEuler[7][0] =   -90.000000*deg;   // 24    South
  ModuleEuler[7][1] =   121.717475*deg;
  ModuleEuler[7][2] =  -126.000000*deg;
  ModuleEuler[8][0] =   180.000000*deg;   // 27    South
  ModuleEuler[8][1] =   148.282526*deg;
  ModuleEuler[8][2] =   126.000000*deg;
  ModuleEuler[9][0] =   180.000000*deg;   // 28    South
  ModuleEuler[9][1] =   148.282526*deg;
  ModuleEuler[9][2] =  -162.000000*deg;
  ModuleEuler[10][0] =  -90.000000*deg;   // 22    Split
  ModuleEuler[10][1] =  121.717475*deg;
  ModuleEuler[10][2] =   90.000000*deg;
  ModuleEuler[11][0] =  180.000000*deg;   // 29    Split
  ModuleEuler[11][1] =  148.282526*deg;
  ModuleEuler[11][2] =  -90.000000*deg;
  ModuleEuler[12][0] =   90.000000*deg;   //  5    North
  ModuleEuler[12][1] =   58.282526*deg;
  ModuleEuler[12][2] =  -18.000000*deg;
  ModuleEuler[13][0] =   90.000000*deg;   //  6    North
  ModuleEuler[13][1] =   58.282526*deg;
  ModuleEuler[13][2] =   54.000000*deg;
  ModuleEuler[14][0] =  -31.750000*deg;   // 10    North
  ModuleEuler[14][1] =   90.000000*deg;
  ModuleEuler[14][2] =  -72.000000*deg;
  ModuleEuler[15][0] =   31.750000*deg;   // 11    North
  ModuleEuler[15][1] =   90.000000*deg;
  ModuleEuler[15][2] =  -36.000000*deg;
  ModuleEuler[16][0] =   31.750000*deg;   // 13    North
  ModuleEuler[16][1] =   90.000000*deg;
  ModuleEuler[16][2] =   36.000000*deg;
  ModuleEuler[17][0] =  -31.750000*deg;   // 14    North
  ModuleEuler[17][1] =   90.000000*deg;
  ModuleEuler[17][2] =   72.000000*deg;
  ModuleEuler[18][0] =  -90.000000*deg;   // 20    North
  ModuleEuler[18][1] =  121.717475*deg;
  ModuleEuler[18][2] =  -54.000000*deg;
  ModuleEuler[19][0] =  -90.000000*deg;   // 21    North
  ModuleEuler[19][1] =  121.717475*deg;
  ModuleEuler[19][2] =   18.000000*deg;
  ModuleEuler[20][0] =  180.000000*deg;   // 25    North
  ModuleEuler[20][1] =  148.282526*deg;
  ModuleEuler[20][2] =  -18.000000*deg;
  ModuleEuler[21][0] =  180.000000*deg;   // 26    North
  ModuleEuler[21][1] =  148.282526*deg;
  ModuleEuler[21][2] =   54.000000*deg;

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

  PosNotch.setX(0.*cm);
  PosNotch.setY(-13.287*2.54*cm);
  PosNotch.setZ(14.646*2.54*cm+notchThickness/2.0);

  PosPlane.setX(0.*cm);
  PosPlane.setY(0.*cm);
  PosPlane.setZ(-Rmax+19.567*2.54*cm);

}

Gretina_NSCL_Shell::~Gretina_NSCL_Shell()
{}

G4int Gretina_NSCL_Shell::FindMaterials()
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


void Gretina_NSCL_Shell::Placement(G4String status)
{
  if ( status != "full" && status != "north" && status != "south" ){
    G4cout << "Shell status " << status << " is not defined." << G4endl;
    return;
  }

  if( FindMaterials() ) return;

  G4RunManager* runManager               = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector = (DetectorConstruction*) runManager->GetUserDetectorConstruction();
  
  //////////////////////////////////////////////////
  // The GRETINA mounting shell
  //////////////////////////////////////////////////

  G4double Phi0 =   0.*deg;
  G4double dPhi = 360.*deg;
  if(status == "north"){
    Phi0 = -90.*deg;
    dPhi = 180.*deg;
  } else if(status == "south"){
    Phi0 =  90.*deg;
    dPhi = 180.*deg;
  }
  G4Sphere* solidShell = new G4Sphere( "solidShell", Rmin, Rmax, Phi0, dPhi, 0., 180.*deg);

  // Beam Port  
  G4Tubs* beamPort = new G4Tubs( "beamPort", 0., smallPortRadius, 1.1*Rmax, 0.*deg, 360.*deg);

  G4SubtractionSolid* shell = new G4SubtractionSolid("Shell", solidShell, beamPort, G4Transform3D(Rot,G4ThreeVector(0.0,0.0,0.0)));

  // Mounting-arm Ports
  G4Tubs* mountPort = new G4Tubs( "mountPort", 0., mountPortRadius, 1.1*Rmax, 0.*deg, 360.*deg);

  Rot.rotateY(90.0*deg);     

  shell = new G4SubtractionSolid ("Shell",shell,mountPort,G4Transform3D(Rot,G4ThreeVector(0.0,0.0,0.0)));

  // S800 Port
  G4Tubs* S800Port = new G4Tubs( "S800Port", 0., S800PortRadius, 4*(Rmax-Rmin)/2., 0., 360.*deg);

  Rot.rotateY(90.0*deg);

  shell = new G4SubtractionSolid ("Shell",shell,S800Port,G4Transform3D(Rot,G4ThreeVector(0.0,0.0,(Rmax+Rmin)/2))); 

  //  G4Tubs* Dewar = new G4Tubs("Dewar", 0, .8*modulePortRadius, .7*(Rmin)/2., 0.*deg, 360.*deg);

  //  G4LogicalVolume* logicDewar = new G4LogicalVolume(Dewar, matShell, "Dewar_log", 0, 0, 0 );

  // Subtract the small ports from the shell.
  G4Tubs* smallPort = new G4Tubs("smallPort", 0, smallPortRadius, 1.5*(Rmax-Rmin)/2., 0., 360.*deg);

  G4int iMin = 0;
  G4int iMax = 10;
  if(status == "north")
    iMin = 5;
  if(status == "south")
    iMax = 5;
  for(G4int i = iMin; i<iMax; i++){
    Rot = G4RotationMatrix::IDENTITY;
    Rot.rotateY( PosSP[i].getTheta() );
    Rot.rotateZ( PosSP[i].getPhi() );
    shell = new G4SubtractionSolid ("Shell",shell,smallPort,G4Transform3D(Rot,PosSP[i]));
  }

  // Subtract the module ports from the shell.
  G4Tubs* modulePort = new G4Tubs("modulePort", 0, modulePortRadius, 1.5*(Rmax-Rmin)/2., 0., 360.*deg);

  iMin = 0;
  iMax = 22;
  if(status == "north")
    iMin = 10;
  if(status == "south")
    iMax = 12;
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


  //  new G4PVPlacement(G4Transform3D(Rot,Pos*.86), logicDewar, "LNDewar", theDetector->HallLog(), false, 0 );

  //  new G4PVPlacement(G4Transform3D(Rot,Pos*.86), logicDewar, "LNDewar", theDetector->HallLog(), false, 0 );

  //  new G4PVPlacement(G4Transform3D(Rot,Pos*.86), logicDewar, "LNDewar", theDetector->HallLog(), false, 0 );

  //  new G4PVPlacement(G4Transform3D(Rot,Pos*.86), logicDewar, "LNDewar", theDetector->HallLog(), false, 0 );

  //  new G4PVPlacement(G4Transform3D(Rot,Pos*.86), logicDewar, "LNDewar", theDetector->HallLog(), false, 0 );

  //  new G4PVPlacement(G4Transform3D(Rot,Pos*.86), logicDewar, "LNDewar", theDetector->HallLog(), false, 0 );

  //  new G4PVPlacement(G4Transform3D(Rot,Pos*.86), logicDewar, "LNDewar", theDetector->HallLog(), false, 0 );

  G4Box* s800notch = new G4Box("s800notch", notchWidth/2.0, 40*cm, notchThickness/2.0);

  Rot = G4RotationMatrix::IDENTITY;

  shell = new G4SubtractionSolid ("Shell",shell,s800notch,G4Transform3D(Rot,PosNotch));

  G4Tubs* s800plane = new G4Tubs( "s800plane", 0., 1.1*Rmax, Rmax, 0., 360.*deg);

  Rot = G4RotationMatrix::IDENTITY;

  //  G4LogicalVolume* logicPlane = new G4LogicalVolume(s800plane, matShell, "Plane_log", 0, 0, 0 );

  //  new G4PVPlacement(G4Transform3D(Rot,PosPlane), logicPlane, "Plane_log", theDetector->HallLog(), false, 0 );

  //  G4IntersectionSolid* Shell = new G4IntersectionSolid ("Shell",shell,s800plane,G4Transform3D(Rot,PosPlane));

  G4LogicalVolume* logicShell = new G4LogicalVolume(shell, matShell, "Shell_log", 0, 0, 0 );

  new G4PVPlacement(0, G4ThreeVector(), "mountingShell", logicShell, theDetector->HallPhys(), false, 0 );

  G4VisAttributes *pVA = new G4VisAttributes( G4Colour(0.0, 1.0, 1.0) );
  logicShell->SetVisAttributes( pVA );

  G4cout << "Constructed the " << status << " shell." << G4endl;
  G4cout << "  Shell radius: " << Rmin << " -- " << Rmax << G4endl;
  G4cout << "  Shell material: " << matShell->GetName() << G4endl;

}
