#ifndef SCANNING
#include "FDS.hh"

FDS::FDS()
{
  myMessenger = new FDS_Messenger(this);
  cloverEulerFile = "";
  shieldEulerFile = "";
  labrEulerFile   = "";
  cloverOuterDL   = 0.5*mm;
  cloverCoaxialDL = 0.5*mm;
}

FDS::~FDS()
{
  delete  myMessenger;
}


void FDS::Placement(G4LogicalVolume* experimentalHall_log)
{
  if(cloverEulerFile != "")
    ReadCloverEulerFile();
  if(shieldEulerFile != "")
    ReadShieldEulerFile();
  if(labrEulerFile != "")
    ReadLaBrEulerFile();
  
  // The detector placement below follows Gretina_Array::PlaceTheClusters()
  
  G4RunManager* runManager = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector
    = (DetectorConstruction*) runManager->GetUserDetectorConstruction();
  
  // Clover Detectors

  G4cout << G4endl << "Placing clovers ... " << G4endl;
  
  G4int nCl, nEa;
  CeulerAngles* pEc = NULL;

  for(nEa = 0; nEa < nCloverEuler; nEa++) {
    pEc = &cloverEuler[nEa];
    nCl = pEc->whichGe;  // (cluster type: N/A here)
    //    G4cout << "##### nCl = pEc->whichGe = " << nCl << G4endl;
    if(nCl < 0) continue;

    materials = new Materials();
    Clover_Detector* aClover
      = new Clover_Detector(experimentalHall_log, materials, "FDS");
    aClover->setTheta(pEc->th);
    aClover->setPhi(pEc->ph);
    aClover->setPsi(pEc->ps);
    aClover->setX(pEc->trasl.x());
    aClover->setY(pEc->trasl.y());
    aClover->setZ(pEc->trasl.z());
    aClover->setCode(31000 + 4*(pEc->numPhys) - 1);
    aClover->setOuterDLThickness(cloverOuterDL);
    aClover->setCoaxialDLThickness(cloverCoaxialDL);
    aClover->Construct();
    aClover->MakeSensitive( theDetector->GetGammaSD() );

    G4cout << " Placed clover " << pEc->numPhys << ", DetCode = "
	   << 31000 + 4*(pEc->numPhys) - 1 << G4endl;
    
  }
  
  // Clover Shields   (... coming real soon now)

  // LaBr Detectors   (... coming real soon now)

}

void FDS::ReadCloverEulerFile()
{
  // Following the euler file approach used in the Gretina_Array class.
  
  FILE      *fp;
  char      line[256];
  G4int     lline, i1, i2;
  float     pps, th, ph, x, y, z;
  
  if( (fp = fopen(cloverEulerFile, "r")) == NULL) {
    G4cout << "\nError opening Euler-angle file " << cloverEulerFile
	   << G4endl;
    exit(EXIT_FAILURE);
  }

  cloverEuler.clear();

  G4cout << "\nReading Euler angles for the clovers from file "
	 << cloverEulerFile << " ..." << G4endl;
  nCloverEuler = 0;
  
  G4RotationMatrix  rm;
  CeulerAngles     *pEa = NULL;

  while(fgets(line, 255, fp) != NULL) {
    lline = strlen(line);
    if(lline < 2) continue;
    if(line[0] == '#') continue;
    if(sscanf(line,"%d %d %f %f %f %f %f %f", &i1, &i2, &pps, &th, &ph, &x, &y, &z) != 8)
      break;
    cloverEuler.push_back( CeulerAngles() );
    pEa = &cloverEuler[nCloverEuler];
    
    pEa->numPhys = i1;
    pEa->whichGe = i2;

    pEa->rotMat.set( 0, 0, 0 );
    pEa->rotMat.rotateZ(((G4double)pps)*deg);
    pEa->rotMat.rotateY(((G4double)th)*deg);
    pEa->rotMat.rotateZ(((G4double)ph)*deg);
    
    pEa->ps      = ((G4double)pps)*deg;
    pEa->th      = ((G4double)th)*deg;
    pEa->ph      = ((G4double)ph)*deg;
    
    pEa->trasl   = G4ThreeVector( ((G4double)x), ((G4double)y), ((G4double)z) ) * mm;
    
    pEa->pTransf = new G4Transform3D( pEa->rotMat, pEa->trasl );
    
    nCloverEuler++;
  }

  fclose(fp);
  G4cout << nCloverEuler << " Euler angles read." << G4endl;
}

void FDS::ReadShieldEulerFile()
{
  FILE      *fp;
  char      line[256];
  G4int     lline, i1, i2;
  float     pps, th, ph, x, y, z;
  
  if( (fp = fopen(shieldEulerFile, "r")) == NULL) {
    G4cout << "\nError opening Euler-angle file " << shieldEulerFile
	   << G4endl;
    exit(EXIT_FAILURE);
  }

  shieldEuler.clear();

  G4cout << "\nReading Euler angles for the clover shields from file "
	 << shieldEulerFile << " ..." << G4endl;
  nShieldEuler = 0;
  
  G4RotationMatrix  rm;
  CeulerAngles     *pEa = NULL;

  while(fgets(line, 255, fp) != NULL) {
    lline = strlen(line);
    if(lline < 2) continue;
    if(line[0] == '#') continue;
    if(sscanf(line,"%d %d %f %f %f %f %f %f", &i1, &i2, &pps, &th, &ph, &x, &y, &z) != 8)
      break;
    shieldEuler.push_back( CeulerAngles() );
    pEa = &shieldEuler[nShieldEuler];
    
    pEa->numPhys = i1;
    pEa->whichGe = i2;

    pEa->rotMat.set( 0, 0, 0 );
    pEa->rotMat.rotateZ(((G4double)pps)*deg);
    pEa->rotMat.rotateY(((G4double)th)*deg);
    pEa->rotMat.rotateZ(((G4double)ph)*deg);
    
    pEa->ps      = ((G4double)pps)*deg;
    pEa->th      = ((G4double)th)*deg;
    pEa->ph      = ((G4double)ph)*deg;
    
    pEa->trasl   = G4ThreeVector( ((G4double)x), ((G4double)y), ((G4double)z) ) * mm;
    
    pEa->pTransf = new G4Transform3D( pEa->rotMat, pEa->trasl );
    
    nShieldEuler++;
  }

  fclose(fp);
  G4cout << nShieldEuler << " Euler angles read." << G4endl;
}

void FDS::ReadLaBrEulerFile()
{
  FILE      *fp;
  char      line[256];
  G4int     lline, i1, i2;
  float     pps, th, ph, x, y, z;
  
  if( (fp = fopen(labrEulerFile, "r")) == NULL) {
    G4cout << "\nError opening Euler-angle file " << labrEulerFile
	   << G4endl;
    exit(EXIT_FAILURE);
  }

  labrEuler.clear();

  G4cout << "\nReading Euler angles for the clover shields from file "
	 << labrEulerFile << " ..." << G4endl;
  nLaBrEuler = 0;
  
  G4RotationMatrix  rm;
  CeulerAngles     *pEa = NULL;

  while(fgets(line, 255, fp) != NULL) {
    lline = strlen(line);
    if(lline < 2) continue;
    if(line[0] == '#') continue;
    if(sscanf(line,"%d %d %f %f %f %f %f %f", &i1, &i2, &pps, &th, &ph, &x, &y, &z) != 8)
      break;
    labrEuler.push_back( CeulerAngles() );
    pEa = &labrEuler[nLaBrEuler];
    
    pEa->numPhys = i1;
    pEa->whichGe = i2;

    pEa->rotMat.set( 0, 0, 0 );
    pEa->rotMat.rotateZ(((G4double)pps)*deg);
    pEa->rotMat.rotateY(((G4double)th)*deg);
    pEa->rotMat.rotateZ(((G4double)ph)*deg);
    
    pEa->ps      = ((G4double)pps)*deg;
    pEa->th      = ((G4double)th)*deg;
    pEa->ph      = ((G4double)ph)*deg;
    
    pEa->trasl   = G4ThreeVector( ((G4double)x), ((G4double)y), ((G4double)z) ) * mm;
    
    pEa->pTransf = new G4Transform3D( pEa->rotMat, pEa->trasl );
    
    nLaBrEuler++;
  }

  fclose(fp);
  G4cout << nLaBrEuler << " Euler angles read." << G4endl;
}

#endif
