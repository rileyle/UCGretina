#ifndef SCANNING
#include "FDS.hh"

FDS::FDS()
{
  myMessenger = new FDS_Messenger(this);
  fdsCloverEulerFile = "";
  fdsShieldEulerFile = "";
  fdsLaBrEulerFile   = "";
}

FDS::~FDS()
{
  delete  myMessenger;
}


void FDS::Placement(G4LogicalVolume* experimentalHall_log)
{
  if(fdsCloverEulerFile != "")
    ReadFDSCloverEulerFile();
  if(fdsShieldEulerFile != "")
    ReadFDSShieldEulerFile();
  if(fdsLaBrEulerFile != "")
    ReadFDSLaBrEulerFile();
  
  // The detector placement below follows Gretina_Array::PlaceTheClusters()
  
  G4RunManager* runManager = G4RunManager::GetRunManager();
  DetectorConstruction* theDetector
    = (DetectorConstruction*) runManager->GetUserDetectorConstruction();
  
  // Clover Detectors

  G4cout << G4endl << "Placing FDS clovers ... " << G4endl;
  
  //  G4int    nGe, nCl, nEa, nCa, nPg, nSol, nPt, indexP;
  G4int nGe, nCl, nEa;
  CeulerAngles* pEc = NULL;

  for(nEa = 0; nEa < nCloverEuler; nEa++) {
    pEc = &cloverEuler[nEa];
    nCl = pEc->whichGe;  // (cluster type: N/A here)
    //    G4cout << "##### nCl = pEc->whichGe = " << nCl << G4endl;
    if(nCl < 0) continue;
    
    nGe = pEc->numPhys; // (slot number)
    G4cout << "##### slot number = nGe = pEc->numPhys = " << nGe << G4endl;

    materials = new Materials();
    Clover_Detector* aClover
      = new Clover_Detector(experimentalHall_log, materials, "FDS");
    aClover->setTheta(pEc->th);
    aClover->setPhi(pEc->ph);
    aClover->setPsi(pEc->ps);
    aClover->setX(pEc->trasl.x());
    aClover->setY(pEc->trasl.y());
    aClover->setZ(pEc->trasl.z());
    aClover->Construct();
    aClover->MakeSensitive( theDetector->GetGammaSD() );
  }
  
  // Clover Shields   (... coming real soon now)

  // LaBr Detectors   (... coming real soon now)

}

void FDS::ReadFDSCloverEulerFile()
{
  // Following the euler file approach used in the Gretina_Array class.
  
  FILE      *fp;
  char      line[256];
  G4int     lline, i1, i2;
  float     pps, th, ph, x, y, z;
  
  if( (fp = fopen(fdsCloverEulerFile, "r")) == NULL) {
    G4cout << "\nError opening Euler-angle file " << fdsCloverEulerFile
	   << G4endl;
    exit(EXIT_FAILURE);
  }

  cloverEuler.clear();

  G4cout << "\nReading Euler angles for the FDS clovers from file "
	 << fdsCloverEulerFile << " ..." << G4endl;
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

void FDS::ReadFDSShieldEulerFile()
{
  FILE      *fp;
  char      line[256];
  G4int     lline, i1, i2;
  float     pps, th, ph, x, y, z;
  
  if( (fp = fopen(fdsShieldEulerFile, "r")) == NULL) {
    G4cout << "\nError opening Euler-angle file " << fdsShieldEulerFile
	   << G4endl;
    exit(EXIT_FAILURE);
  }

  shieldEuler.clear();

  G4cout << "\nReading Euler angles for the FDS clover shields from file "
	 << fdsShieldEulerFile << " ..." << G4endl;
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

void FDS::ReadFDSLaBrEulerFile()
{
  FILE      *fp;
  char      line[256];
  G4int     lline, i1, i2;
  float     pps, th, ph, x, y, z;
  
  if( (fp = fopen(fdsLaBrEulerFile, "r")) == NULL) {
    G4cout << "\nError opening Euler-angle file " << fdsLaBrEulerFile
	   << G4endl;
    exit(EXIT_FAILURE);
  }

  labrEuler.clear();

  G4cout << "\nReading Euler angles for the FDS clover shields from file "
	 << fdsLaBrEulerFile << " ..." << G4endl;
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
