
#include "Incoming_Beam.hh"


Incoming_Beam::Incoming_Beam()
{
  A=20;
  Z=12;
  Ex=0.;
  //  KEu=100*MeV;
  KEu=0;
  KE=KEu*A;
  dtaFileName = "";
  Ndta = 0;
  dtaMin = 0.;
  dtaMax = 0.;
  dtaBin = 0.;
  Dpp=0.0;
  fcX=0.;
  fcDX=0.;
  fcY=0.;
  fcDY=0.;
  fcZ=-50.*cm;
  maxAta=0*mrad;
  maxBta=0*mrad;
  ata0=0.*mrad;
  bta0=0.*mrad;
}

Incoming_Beam::~Incoming_Beam()
{;}
//---------------------------------------------------------
void Incoming_Beam::Report()
{

  G4cout<<"----> Z of the incoming beam set to  "<<Z<< G4endl;
  G4cout<<"----> A of the incoming beam set to  "<<A<< G4endl;
  G4cout<<"----> Kin. En. of the incoming beam set to "<<
  G4BestUnit(KE,"Energy")<<G4endl;
  if(dtaFileName != ""){
    G4cout<<"----> Relative KE distribution of the incoming beam read from "
	  << dtaFileName << G4endl;
    G4cout<<"----> KE per nucleon corresponding to zero relative KE set to "
	  << G4BestUnit(KEu,"Energy")<<G4endl;
  } else
    G4cout<<"----> KE per nucleon of the incoming beam set to "<<
      G4BestUnit(KEu,"Energy")<<G4endl;
  G4cout<<"----> momentum acceptance for the incoming beam set to  "<<Dpp<< G4endl;
  G4cout<<"----> focal point X position for the incoming beam set to  "<<G4BestUnit(fcX,"Length")<< G4endl;
  G4cout<<"----> focal point DX size for the incoming beam set to  "<<G4BestUnit(fcDX,"Length")<< G4endl;
  G4cout<<"----> focal point Y position for the incoming beam set to  "<<G4BestUnit(fcY,"Length")<< G4endl;
  G4cout<<"----> focal point DY size for the incoming beam set to  "<<G4BestUnit(fcDY,"Length")<< G4endl;
  G4cout<<"----> focal point Z position for the incoming beam set to  "<<G4BestUnit(fcZ,"Length")<< G4endl;
  G4cout<<"----> dispersive direction angular divergence for the incoming beam set to  "<<maxAta/mrad<<" mrad = "<<maxAta/deg<<" deg"<< G4endl;
  G4cout<<"----> non dispersive direction angular divergence for the incoming beam set to  "<<maxBta/mrad<<" mrad = "<<maxBta/deg<<" deg"<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setA(G4int Ain)
{

  A=Ain;
  //  G4cout<<"----> A of the incoming beam set to  "<<A<< G4endl;
  
}
//---------------------------------------------------------
void Incoming_Beam::setZ(G4int Zin)
{

  Z=Zin;
  //  G4cout<<"----> Z of the incoming beam set to  "<<Z<< G4endl;
  
}
//---------------------------------------------------------
void Incoming_Beam::setDpp(G4double d)
{

  Dpp=d;
  //  G4cout<<"----> momentum acceptance for the incoming beam set to  "<<Dpp<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setfcX(G4double d)
{

  fcX=d;
  //  G4cout<<"----> focal point X position for the incoming beam set to  "<<G4BestUnit(fcX,"Length")<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setfcDX(G4double d)
{

  fcDX=d;
  //  G4cout<<"----> focal point DX size for the incoming beam set to  "<<G4BestUnit(fcDX,"Length")<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setfcDY(G4double d)
{

  fcDY=d;
  //  G4cout<<"----> focal point DY size for the incoming beam set to  "<<G4BestUnit(fcDY,"Length")<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setfcY(G4double d)
{

  fcY=d;
  //  G4cout<<"----> focal point Y position for the incoming beam set to  "<<G4BestUnit(fcY,"Length")<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setfcZ(G4double d)
{

  fcZ=d;
  //  G4cout<<"----> focal point Z position for the incoming beam set to  "<<G4BestUnit(fcZ,"Length")<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setmaxAta(G4double d)
{

  maxAta=d;
  //  G4cout<<"----> dispersive direction angular divergence for the incoming beam set to  "<<maxAta/mrad<<" mrad = "<<maxAta/deg<<" deg"<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setmaxBta(G4double d)
{

  maxBta=d;
  //  G4cout<<"----> non dispersive direction angular divergence for the incoming beam set to  "<<maxBta/mrad<<" mrad = "<<maxBta/deg<<" deg"<< G4endl;
}
//---------------------------------------------------------
void Incoming_Beam::setKE(G4double KEin)
{

  KE=KEin;
  KEu=KE/A;
  //G4cout<<"----> Kin. En. of the incoming beam set to "<<
  // G4BestUnit(KE,"Energy")<<G4endl;
  //G4cout<<"----> Kin. En. per nucleon of the incoming beam set to "<<
  // G4BestUnit(KEu,"Energy")<<G4endl;
}
//-----------------------------------------------------------------
void Incoming_Beam::setKEu(G4double KEuin)
{

  KEu=KEuin;
  KE=KEu*A;
  //G4cout<<"----> Kin. En. of the incoming beam set to "<<
  // G4BestUnit(KE,"Energy")<<G4endl;
  //G4cout<<"----> Kin. En. per nucleon of the incoming beam set to "<<
  // G4BestUnit(KEu,"Energy")<<G4endl; 
}
//-----------------------------------------------------------------
void Incoming_Beam::setDTAFile(G4String fileName)
{
  char line[1000];

  dtaFileName = fileName;
  std::ifstream dtaFile;
  dtaFile.open(dtaFileName);
  if (!dtaFile.is_open()) {
    G4cout << "Unable to open " << dtaFileName << G4endl;
    exit(EXIT_FAILURE);
  }
 
  dtaFile >> dtaMin >> dtaMax >> dtaBin;
  dtaFile.getline(line,1000);  // Advance to next line.

  G4double d;
  Ndta = 0;
  while(dtaFile >> d){
    dtaFile.getline(line,1000);  // Advance to next line.
    dta[Ndta] = d;
    Ndta++;
  }

  dtaFile.close();

}
//---------------------------------------------------------
G4ThreeVector Incoming_Beam::getDirection()
{
  G4ThreeVector direction;
  G4double x,y,z,a,b,r,phi;
  //  r1=G4UniformRand();
  //   r2=G4UniformRand();
  
  //   z=r1+cos(maxTh)*(1-r1);
  //   x=sqrt(1-z*z)*sin(2*pi*r2);
  //   y=sqrt(1-z*z)*cos(2*pi*r2);

  phi=G4UniformRand()*8.*atan(1.);
  r=G4UniformRand()+G4UniformRand();
  if(r>=1) r=-(r-2.);

  a=r*cos(phi)*maxAta;
  b=r*sin(phi)*maxBta; 
  z=1./sqrt(1.+tan(a)*tan(a)+tan(b)*tan(b));
  y=z*tan(b);
  x=z*tan(a);
  direction.setX(x);
  direction.setY(y);
  direction.setZ(z);
  direction.rotateY(ata0);
  direction.rotateX(-bta0);
  return direction;
  
}
//---------------------------------------------------------
G4ThreeVector Incoming_Beam::getPosition()
{
  G4ThreeVector position;
  G4double x,y;
  G4double r,phi;
  
  phi=G4UniformRand()*8.*atan(1.);
  r=G4UniformRand()+G4UniformRand();
  if(r>=1) r=-(r-2.);

  x=fcX+r*cos(phi)*fcDX/2.;
  y=fcY+r*sin(phi)*fcDY/2.;

  //At emission point!!! (Macro file command names are misleading.)
  position.setX(x);
  position.setY(y);
  position.setZ(fcZ);
  return position;
  
}

//---------------------------------------------------------
G4double Incoming_Beam::getKE(G4ParticleDefinition *ion)
{
  G4DynamicParticle dynamic;
  G4ThreeVector momentum_vector;
  G4double momentum;
  G4double rand;
  G4double ke, KinEne;
  //  momentum_vector=0;    //LR (Change to CLHEP library Hep3Vector)
  momentum_vector.setX(0.); //LR
  momentum_vector.setY(0.); //LR
  momentum_vector.setZ(1.);

  if(Ndta > 0){
    CLHEP::RandGeneral randDTA( dta, Ndta );
    ke = KE*(1 + dtaMin/100. + randDTA.fire()*(dtaMax-dtaMin)/100.);
  } else {
    ke = KE;
  }

  dynamic=G4DynamicParticle(ion,momentum_vector,ke);
  momentum=dynamic.GetTotalMomentum();
  rand=G4UniformRand()-0.5;
  momentum*=(1+rand*Dpp);
  momentum_vector.setMag(momentum);
  //  dynamic.SetMomentum(momentum); //LR (Change to CLHEP library Hep3Vector)
  dynamic.SetMomentum(momentum_vector);
  KinEne=dynamic.GetKineticEnergy();

  return KinEne;
  
}
