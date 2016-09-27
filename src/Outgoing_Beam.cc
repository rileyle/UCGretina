#include "Outgoing_Beam.hh"
#include "G4DecayTable.hh"
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4ProcessManager.hh"

Outgoing_Beam::Outgoing_Beam()
{
  DZ=0;
  DA=0;
  Ex=1.0*MeV;
  lvlDataFileName = "";
  xsectFileName = "";
  TarEx=0.*keV;
  TarA = 1;
  TarZ = 1;
  TFrac=0.;

  targetExcitation=false;
  reacted=false;
  source=false;
  sigma_a=0.;
  sigma_b=0.;
  theta_max=180.*deg;
  theta_min=0.;
  theta_bin=0.;
  twopi=8.*atan(1.);
  beamIn = NULL;
}

Outgoing_Beam::~Outgoing_Beam()
{
  ;
}

void Outgoing_Beam::defaultIncomingIon(Incoming_Beam *bi)
{
  if (bi == NULL) {
    return;
  }
  beamIn = bi;
}

void Outgoing_Beam::setDecayProperties()
{
  if (beamIn == NULL) {
    G4cerr << "Can not set decay properties as incoming beam has not been set" << G4endl;
    exit(EXIT_FAILURE);
  }

  Zin = beamIn->getZ();
  Ain = beamIn->getA();

  // Load the particle table with ground states
  beam  = G4IonTable::GetIonTable()->GetIon(Zin,Ain,0.);
  ionGS = G4IonTable::GetIonTable()->GetIon(Zin+DZ,Ain+DA,0.);

  if(TarA == 1 && TarZ ==1)
    tarIn = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  else if (TarA == 2 && TarZ ==1)
    tarIn = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");
  else if (TarA == 3 && TarZ ==1)
    tarIn = G4ParticleTable::GetParticleTable()->FindParticle("triton");
  else if (TarA == 3 && TarZ ==2)
    tarIn = G4ParticleTable::GetParticleTable()->FindParticle("He3");
  else if (TarA == 1 && TarZ ==0)
    tarIn = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  else
    tarIn = G4IonTable::GetIonTable()->GetIon(TarZ,    TarA,    0.);

  if(DA == 0){
    tarOut   = tarIn;
    tarOutGS = tarIn;
  } else if(TarA-DA == 1 && TarZ-DZ ==1){
    tarOut = G4ParticleTable::GetParticleTable()->FindParticle("proton");
    tarOutGS = tarOut;
  } else if(TarA-DA == 2 && TarZ-DZ ==1){
    tarOut = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");
    tarOutGS = tarOut;
  } else if(TarA-DA == 3 && TarZ-DZ ==1){
    tarOut = G4ParticleTable::GetParticleTable()->FindParticle("triton");
    tarOutGS = tarOut;
  } else if(TarA-DA == 3 && TarZ-DZ ==2) {
    tarOut = G4ParticleTable::GetParticleTable()->FindParticle("He3");
    tarOutGS = tarOut;
  } else if(TarA-DA == 1 && TarZ-DZ ==0) {
    tarOut = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
    tarOutGS = tarOut;
  } else {
    tarOut = G4IonTable::GetIonTable()->GetIon(TarZ-DZ, TarA-DA, Ex);
    tarOutGS = G4IonTable::GetIonTable()->GetIon(TarZ-DZ, TarA-DA, 0.);
  }
  if (ionGS == NULL) {
    G4cerr << "Error: no outgoing ion ground state in particle table "
	   << Zin + DZ << " " << Ain+DA << G4endl;
    exit(EXIT_FAILURE);
  }
  if (tarIn == NULL) {
    if(targetExcitation){
      G4cerr << "Error: no target nucleus in particle table "
	     << TarZ << " " << TarA << G4endl;
      exit(EXIT_FAILURE);
    } else {
      G4cerr << "Warning: no target nucleus in particle table "
	     << TarZ << " " << TarA << G4endl;
    }
  }
  if (tarOut == NULL) {
    if(targetExcitation){
      G4cerr << "Error: no target-like product in particle table "
	     << TarZ-DZ << " " << TarA-DA << G4endl;
      exit(EXIT_FAILURE);
    } else {
      G4cerr << "Warning: no target-like product in particle table "
	     << TarZ-DZ << " " << TarA-DA << G4endl;
    }
  }
  if (tarOutGS == NULL) {
    if(targetExcitation){
      G4cerr << "Error: no target-like ground-state product in particle table "
	     << TarZ-DZ << " " << TarA-DA << G4endl;
      exit(EXIT_FAILURE);
    } else {
      G4cerr << "Warning: no target-like ground-state product in particle table "
	     << TarZ-DZ << " " << TarA-DA << G4endl;
    }
  }

  m1 = beam->GetPDGMass();
  m2 = tarIn->GetPDGMass();
  if(targetExcitation){  
    m3 = ionGS->GetPDGMass();
    m4 = tarOut->GetPDGMass();
  } else {
    m3 = ionGS->GetPDGMass();
    if (tarOutGS == NULL)
      m4 = 0.;
    else
      m4 = tarOutGS->GetPDGMass();
  }

  // Load the particle table with excited states
  G4int Z,A;

  if(targetExcitation){
    Z = tarOut->GetAtomicNumber();
    A = tarOut->GetAtomicMass();
  } else {
    Z = ionGS->GetAtomicNumber();
    A = ionGS->GetAtomicMass();
  }

  if(lvlDataFileName == ""){
    G4cerr << "Error: level data file name not set." << G4endl;
    exit(EXIT_FAILURE);
  }
  
  G4NuclearLevelData* levelData = G4NuclearLevelData::GetInstance();
  levelData->AddPrivateData(Z, A, lvlDataFileName);
  const G4LevelManager* levelManager = levelData->GetLevelManager(Z, A);
  G4int Nentries = levelManager->NumberOfTransitions()+1;
  for(G4int i = 1; i < Nentries; i++){ // Excited states
    // G4cout << "Level " << i
    // 	   << " energy = " << levelManager->LevelEnergy(i)
    // 	   << G4endl;
    
    G4ParticleDefinition* excitedState
      = G4IonTable::GetIonTable()->GetIon(Z,A,levelManager->LevelEnergy(i));

    excitedState->SetPDGStable(false);
    excitedState->SetPDGLifeTime(levelManager->LifeTime(i));
  }

  //  G4IonTable::GetIonTable()->DumpTable();

}

//---------------------------------------------------------
void Outgoing_Beam::ScanInitialConditions(const G4Track & aTrack)
{
 
  dirIn=aTrack.GetMomentumDirection();
  posIn=aTrack.GetPosition();
  pIn=aTrack.GetMomentum();
  KEIn=aTrack.GetDynamicParticle()->GetKineticEnergy();
  Ain=aTrack.GetDynamicParticle()->GetDefinition()->GetAtomicMass();
  Zin=aTrack.GetDynamicParticle()->GetDefinition()->GetAtomicNumber();
  ReactionFlag=-1;
  if(aTrack.GetVolume()->GetLogicalVolume()->GetName()=="target_log")
    ReactionFlag=0;
  ThresholdFlag = 1;
  if( KEIn < Ex/2/m2*(m1+m2+m3+m4) )
    ThresholdFlag = 0;
  ET = KEIn + m1 + m2;                                      // Total E
  p1 = sqrt( (KEIn + m1)*(KEIn + m1) - m1*m1 );          // Incoming p
  // Lab-frame scattering angle limit
  sin2theta3_max = 
    ( (ET*ET - p1*p1 + m3*m3 - m4*m4)*(ET*ET - p1*p1 + m3*m3 - m4*m4)
      -4*m3*m3*((m1 + m2)*(m1 + m2) + 2*m2*KEIn) )/( 4*m3*m3*p1*p1 );

}

//---------------------------------------------------------
G4ThreeVector Outgoing_Beam::ReactionPosition()
{
 
  posOut=posIn;
  posOut.setZ(posIn.getZ()+eps);
  //  G4cout<<"Reaction took place at "<<G4BestUnit(posOut.getZ(),"Length")<<G4endl;

  return posOut;

}
//---------------------------------------------------------
G4DynamicParticle* Outgoing_Beam::ReactionProduct()
{

  G4int Zout, Aout;
  G4double excitationEnergy = 0.;
  if(targetExcitation){
    Zout = TarZ - DZ;
    Aout = TarA - DA;
    excitationEnergy = TarEx;
  } else {
    Zout = Zin + DZ;
    Aout = Ain + DA;
    excitationEnergy = Ex;
  }

  G4ParticleDefinition* product
    = G4IonTable::GetIonTable()->GetIon(Zout, Aout, excitationEnergy);

  G4DynamicParticle* aReactionProduct
    = new G4DynamicParticle(product, GetOutgoingMomentum());

  return aReactionProduct;
}

//---------------------------------------------------------
G4ThreeVector Outgoing_Beam::GetOutgoingMomentum()
{
  static const G4ThreeVector ez(0.,0.,1.);

  G4ThreeVector ax,ppOut;
  G4double      theta3;

  // Lab-frame scattering angle ================================================

  theta3=GetDTheta();
  // If theta3 is beyond the limit dictated by the kinematics
  // or if an angle cut is specified ...
  while( sin(theta3)*sin(theta3) > sin2theta3_max
	 || (theta3 < theta_min)
	 || (theta3 > theta_max) ){
    theta3=GetDTheta();
  }

  // Relativistic kinematics ===================================================
  //       Baldin et al, Kinematics of Nuclear Reactions, Pergamon (1961)

  G4double E3Lab = 1/(ET*ET - p1*p1*cos(theta3)*cos(theta3))// Beam-like product
    *( ET*( m2*(KEIn + m1) + (m1*m1+m2*m2+m3*m3-m4*m4)/2. )
       + p1*cos(theta3)*sqrt( ( m2*(KEIn+m1) + (m1*m1+m2*m2-m3*m3-m4*m4)/2. )
			     *( m2*(KEIn+m1) + (m1*m1+m2*m2-m3*m3-m4*m4)/2. ) 
			     - m3*m3*m4*m4 
			     - p1*p1*m3*m3*sin(theta3)*sin(theta3) ) );

  G4double p3Lab = sqrt( E3Lab*E3Lab - m3*m3 );

  G4double p4Lab = sqrt( (ET-E3Lab)*(ET-E3Lab) - m4*m4);  // Target-like product

  G4double theta4 = asin( p3Lab/p4Lab * sin(theta3) );

  G4double pLab = p3Lab;       // Construct pLab vector of the Beam-like product
  G4double theta = theta3;
  if( targetExcitation ){    // Construct pLab vector of the Target-like product
    pLab = p4Lab;
    theta = theta4;
  }

  // Set the magnitude of the outgoing momentum ================================

  ppOut = pIn;
  if( ppOut.mag() > 0)
      ppOut.setMag( pLab );

  // Set the direction of the outgoing momentum ================================

  ax=pIn.cross(ez);
  if (pIn.mag2() == 0.) {
    return ppOut;
  }
  ax.rotate(pIn,G4UniformRand()*twopi);
  
  if (ax.mag2() == 0.) {
    return ppOut;
  }
  ppOut.rotate(ax,theta);

  //  G4cout << std::setprecision(3) << std::setw(10)
  //	 << "ppOut.theta() = " << ppOut.theta()/deg 
  //	 << "  ata = " << -asin(ppOut.getY()/ppOut.mag())/mrad
  //	 << "  bta = " << -asin(ppOut.getX()/ppOut.mag())/mrad
  //	 << G4endl;
  
  return ppOut;
}

//---------------------------------------------------------
G4double Outgoing_Beam::GetDTheta()
{

  G4double theta;
  if(xsectFileName == ""){
    G4double theta_a,theta_b;
    theta_a = CLHEP::RandGauss::shoot(0,sigma_a);
    theta_b = CLHEP::RandGauss::shoot(0,sigma_b);
    theta = sqrt(theta_a*theta_a+theta_b*theta_b);
  } else {
    CLHEP::RandGeneral randTheta( Xsect, Nxsect );
    theta = theta_min + randTheta.fire()*(theta_max-theta_min);
  }

  return theta;
}
//---------------------------------------------------------
void Outgoing_Beam::setXsectFile(G4String fileName)
{
 
  char line[1000];

  xsectFileName = fileName;
  std::ifstream xsectFile;
  xsectFile.open(xsectFileName);

  xsectFile >> theta_min >> theta_max >> theta_bin;
  xsectFile.getline(line,1000);  // Advance to next line.

  theta_min *= deg;
  theta_max *= deg;
  theta_bin *= deg;

  G4double x;
  Nxsect = 0;
  while(xsectFile >> x){
    xsectFile.getline(line,1000);  // Advance to next line.
    Xsect[Nxsect] = x*sin( theta_min + (Nxsect+0.5)*theta_bin );
    Nxsect++;
  }

  xsectFile.close();

}
//---------------------------------------------------------
void Outgoing_Beam::Report()
{

  G4cout<<"----> Delta A for the outgoing beam set to  "<<DA<< G4endl;
  G4cout<<"----> Delta Z for the outgoing beam set to  "<<DZ<< G4endl;
  G4cout<<"----> Target A (for kinematics calculations) is set to  "<<TarA<< G4endl;
  G4cout<<"----> Target Z (for kinematics calculations) is set to  "<<TarZ<< G4endl;
  G4cout<<"----> Excitation energy of the outgoing beam set to "<<
    G4BestUnit(Ex,"Energy")<<G4endl;
  G4cout<<"----> Target excitation energy set to "<< G4BestUnit(TarEx,"Energy")<<G4endl;  
  G4cout<<"----> Fraction of target excitations set to "<<TFrac<<G4endl;
  G4cout<<"----> Mass of the incoming beam is " <<beam->GetPDGMass()<<" MeV"<<G4endl;
  G4cout<<"----> Mass of the target is " <<tarIn->GetPDGMass()<<" MeV"<<G4endl;
  G4cout<<"----> Mass of the outgoing beam-like reaction product is " <<ionGS->GetPDGMass()<<" MeV"<<G4endl;
  if(tarOut != NULL)
    G4cout<<"----> Mass of the outgoing target-like reaction product is " <<tarOut->GetPDGMass()<<" MeV"<<G4endl;
  G4cout<<"----> Sigma for ata distribution set to "<<sigma_a<<G4endl;
  G4cout<<"----> Sigma for bta distribution set to "<<sigma_b<<G4endl;

  G4IonTable::GetIonTable()->DumpTable();
}


//---------------------------------------------------------
void Outgoing_Beam::setDA(G4int da)
{

  DA=da;
  // G4cout<<"----> Delta A for the outgoing beam set to  "<<DA<< G4endl; 
}
//---------------------------------------------------------
void Outgoing_Beam::setDZ(G4int dz)
{

  DZ=dz;
  //  G4cout<<"----> Delta Z for the outgoing beam set to  "<<DZ<< G4endl;
  
}
//---------------------------------------------------------
void Outgoing_Beam::setEx(G4double ex)
{

  Ex=ex;
  
  //  G4cout<<"----> Excitation energy of the outgoing beam set to "<<
  //  G4BestUnit(Ex,"Energy")<<G4endl;

}
//---------------------------------------------------------
void Outgoing_Beam::setTarEx(G4double ex)
{

  TarEx=ex;
  
  //G4cout<<"----> Target excitation energy set to "<< G4BestUnit(TarEx,"Energy")<<G4endl;

}
//---------------------------------------------------------
void Outgoing_Beam::setTarA(G4int a)
{

  TarA = a;
  
  //  G4cout<<"----> Target A set to "<< TarA <<G4endl;

}
//---------------------------------------------------------
void Outgoing_Beam::setTarZ(G4int z)
{

  TarZ = z;
  
  //  G4cout<<"----> Target Z set to "<< TarZ <<G4endl;

}
//---------------------------------------------------------
void Outgoing_Beam::setTFrac(G4double ex)
{
  TFrac=ex;
  if(TFrac>1) TFrac=1.;
  if(TFrac<0) TFrac=0.;

  //G4cout<<"----> Fraction of target excitations set to "<<TFrac<<G4endl;

}
