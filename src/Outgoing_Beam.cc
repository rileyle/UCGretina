#include "Outgoing_Beam.hh"
#include "GammaDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4Decay.hh"
#include "G4ProcessManager.hh"

G4Decay Outgoing_Beam::decay;

Outgoing_Beam::Outgoing_Beam()
{
  DZ=0;
  DA=0;
  Ex=1.0*MeV;
  lvlSchemeFileName = "";
  TarEx=547.5*keV;
  TFrac=0.;
  tau=0.0*ns;
  betaDopp=0.33;
  posDopp.setX(0.);
  posDopp.setY(0.);
  posDopp.setZ(0.);

  reacted=false;
  alpha=20.;
  sigma_a=0.;
  sigma_b=0.;
  theta_max=-1.;
  Calc_pmax();
  twopi=8.*atan(1.);
  NQ=1;SQ=0;
  SetUpChargeStates();
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
  ion = G4ParticleTable::GetParticleTable()->GetIon(Zin+DZ,Ain+DA,Ex);
  iongs = G4ParticleTable::GetParticleTable()->GetIon(Zin+DZ,Ain+DA,0.);
  iongs->SetPDGStable(true);
  if (ion == NULL) {
    G4cerr << "Could not find outgoing ion in particle table "
	   << Zin + DZ << " " << Ain+DA << G4endl;
    exit(EXIT_FAILURE);
  }
  if (iongs == NULL) {
    G4cerr << "Could not find outgoing ion in particle table "
	   << Zin + DZ << " " << Ain+DA << G4endl;
    exit(EXIT_FAILURE);
  }

  G4DecayTable *DecTab = NULL;
  GammaDecayChannel *GamDec = NULL;
  G4ProcessManager *pm = NULL;

  if (lvlSchemeFileName == ""){
    G4cout << "Constructing decay properties for Z=" << Zin + DZ
	   << " A=" << Ain + DA << " with excitation " << Ex/keV << " keV" << G4endl;
    G4cout << "Direct gamma decay to the ground state." << G4endl;

    ion->SetPDGStable(false);
    ion->SetPDGLifeTime(tau);

    DecTab = ion->GetDecayTable(); 
    if (DecTab == NULL) {
      DecTab = new G4DecayTable();
      ion->SetDecayTable(DecTab);
    }
    GammaDecayChannel *GamDec = new GammaDecayChannel(-1,ion,1,Ex,0.,theAngularDistribution);
    DecTab->Insert(GamDec);
    // DecTab->DumpInfo();

    // make sure that the ion has the decay process in its manager
    G4ProcessManager *pm = ion->GetProcessManager();
    if (pm == NULL) {
      G4cerr << "Could not find process manager for outgoing ion." << G4endl;
      exit(EXIT_FAILURE);
    }
    pm->AddProcess(&decay,1,-1,4);
    // pm->DumpInfo();

  } else { // Set up intermediate states and decay properties
    G4cout << "Reading level scheme description from " 
	   << lvlSchemeFileName << G4endl;

    G4int nBranch;
    G4double Elevel, meanLife, BR, Exf, a0, a2, a4;
    G4ParticleDefinition* intermediateIon;

    openLvlSchemeFile();

    G4int i = 0;
    while(lvlSchemeFile >> Elevel >> nBranch >> meanLife){
      G4cout << "Constructing decay properties for Z=" << Zin + DZ
	     << " A=" << Ain + DA << " with excitation " << Elevel << " keV" << G4endl;
      for(G4int j = 0; j < nBranch; j++){
	lvlSchemeFile >> BR >> Exf >> a0 >> a2 >> a4;

	if(i == 0)
	  intermediateIon = ion;
	else
	  intermediateIon = G4ParticleTable::GetParticleTable()->GetIon(Zin+DZ,Ain+DA,Elevel*keV);
	if (intermediateIon == NULL) {
	  G4cerr << "Could not find intermediate ion in particle table "
		 << Zin + DZ << " " << Ain+DA << G4endl;
	  exit(EXIT_FAILURE);
	}
	intermediateIon->SetPDGStable(false);
	intermediateIon->SetPDGLifeTime(meanLife*picosecond);

	DecTab = intermediateIon->GetDecayTable(); 
	if (DecTab == NULL) {
	  DecTab = new G4DecayTable();
	  intermediateIon->SetDecayTable(DecTab);
	}

	theAngularDistribution.SetCoeffs(a0,a2,a4);
	theAngularDistribution.Report();

	GamDec = new GammaDecayChannel(-1,intermediateIon,BR,(Elevel-Exf)*keV,Exf*keV,theAngularDistribution);
	DecTab->Insert(GamDec);

	// make sure that the ion has the decay process in its manager
	pm = intermediateIon->GetProcessManager();
	if (pm == NULL) {
	  G4cerr << "Could not find process manager for outgoing ion." << G4endl;
	  exit(EXIT_FAILURE);
	}
	pm->AddProcess(&decay,1,-1,4);

      }
      cout << endl;
      DecTab->DumpInfo();
      //      pm->DumpInfo();

      i++;
    }

    closeLvlSchemeFile();

    //    G4ParticleTable::GetParticleTable()->DumpTable();

  }
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
  tauIn=aTrack.GetProperTime();
  ReactionFlag=-1;
  if(aTrack.GetVolume()->GetLogicalVolume()->GetName()=="target_log")
    ReactionFlag=0;
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
  G4DynamicParticle* aReactionProduct =new G4DynamicParticle(ion,GetOutgoingMomentum());

  //  aReactionProduct->SetProperTime(tauIn);
  // G4cout<<" Proper time set to "<<aReactionProduct->GetProperTime()<<G4endl;

  return aReactionProduct;
}
//---------------------------------------------------------
G4DynamicParticle* Outgoing_Beam::ProjectileGS()
{
  G4DynamicParticle* aReactionProduct =new G4DynamicParticle(iongs,GetOutgoingMomentum());

  //  aReactionProduct->SetProperTime(tauIn);
  // G4cout<<" Proper time set to "<<aReactionProduct->GetProperTime()<<G4endl;

  return aReactionProduct;
}
//---------------------------------------------------------
G4DynamicParticle* Outgoing_Beam::TargetExcitation()
{
  particleTable = G4ParticleTable::GetParticleTable();
  G4DynamicParticle* aReactionProduct =new G4DynamicParticle(particleTable->FindParticle("gamma"),TargetAngularDistribution(),TarEx);

  return aReactionProduct;
}

//---------------------------------------------------------
G4ThreeVector Outgoing_Beam::TargetAngularDistribution()
{
    //TB allow a coulex angular distribution for target excitation (coulex on gold)
    G4ThreeVector direction = G4RandomDirection();
    direction.setTheta(theTargetAngularDistribution.GetRandomAngle());
    //G4cout<<direction.getTheta()<<G4endl;
    return direction;
}
//---------------------------------------------------------
G4ThreeVector Outgoing_Beam::GetOutgoingMomentum()
{
  static const G4ThreeVector ez(0.,0.,1.);

  //? (RFE) differentiate Coulomb scattering vs. Knockout
  G4ThreeVector ax,ppOut;
  G4double      theta;
  // ppin/Ain = ppout/Aout, Aout = Ain + DA
  //ppOut=pIn * (float(Ain) / float(Ain-DA)); // conserve KE/u
  //  ppOut=pIn * (float(Ain) / float(Ain-DA)); // conserve KE/u  RS (Aout=Ain+DA)

  // conserve KE/u LR
  //  G4double tout = KEIn*float(Ain+DA)/float(Ain);
  //  ppOut = pIn;
  //  ppOut.setMag( sqrt(tout*tout + 2*ion->GetPDGMass()*tout) );
  
  // conserve p/u LR
  ppOut=pIn * float(Ain+DA)/float(Ain);

  ax=pIn.cross(ez);
  if (pIn.mag2() == 0.) {
    return ppOut;
  }
  ax.rotate(pIn,G4UniformRand()*twopi);

  theta=GetDTheta();
  //  G4cout << std::setprecision(3) << std::setw(10)
  //	 << "Theta_0 = " << theta/deg << G4endl;
  //If an angle cut is specified ...
  if(theta_max > 0.)
    while(theta > theta_max){
      theta=GetDTheta();
      //      G4cout << std::setprecision(3) << std::setw(10)
      //	     << "***Theta_1 = " << theta/deg << G4endl;

    }

  
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
  //LR G4double delta,d,theta;

  //LR delta=alpha/2.-1./12.;
  //LR d=sqrt(1+4.*delta*pmax*G4UniformRand());
  //LR theta=sqrt(0.5*(-1.+d)/delta);

  G4double theta_a,theta_b,theta;

  theta_a = CLHEP::RandGauss::shoot(0,sigma_a);
  theta_b = CLHEP::RandGauss::shoot(0,sigma_b);
  theta = sqrt(theta_a*theta_a+theta_b*theta_b);

  return theta;
}
//---------------------------------------------------------
void Outgoing_Beam::Calc_pmax()
{

  G4double delta;
  delta=alpha/2.-1./12.;
  pmax=(1.+delta*theta_max*theta_max)*theta_max*theta_max;
    
}
//---------------------------------------------------------
void Outgoing_Beam::Report()
{
  setDecayProperties();

G4cout<<"----> Delta Z for the outgoing beam set to  "<<DZ<< G4endl;
G4cout<<"----> Delta A for the outgoing beam set to  "<<DA<< G4endl;
G4cout<<"----> Excitation energy of the outgoing beam set to "<<
 G4BestUnit(Ex,"Energy")<<G4endl;
G4cout<<"----> Target excitation energy set to "<< G4BestUnit(TarEx,"Energy")<<G4endl;  
G4cout<<"----> Fraction of target excitations set to "<<TFrac<<G4endl;
G4cout<<"----> Lifetime of the excited state for the outgoing beam set to "<<
 G4BestUnit(tau,"Time")<<G4endl; 
//G4cout<<"----> Beta for the outgoing beam set to "<<betaDopp<<G4endl; 
//G4cout<<"----> X position for the outgoing beam Doppler corrections set to "<<posDopp.getX()<<G4endl; 
//G4cout<<"----> Y position for the outgoing beam Doppler corrections set to "<<posDopp.getY()<<G4endl; 
//G4cout<<"----> Z position for the outgoing beam Doppler corrections set to "<<posDopp.getZ()<<G4endl;
//G4cout<<"----> Alpha coefficient for angular distribution of ions Coulomb scattered on the target set to "<<alpha<<G4endl; 
// G4cout<<"----> Maximum scattering angle for ions Coulomb scattered on the target "<<theta_max/deg<<" deg. "<<G4endl; 
 G4cout<<"Sigma for ata............"<<sigma_a<<G4endl;
 G4cout<<"Sigma for bta............"<<sigma_b<<G4endl;
 G4cout<<"----> Number of charge states "<<NQ<<G4endl;
 vector<Charge_State*>::iterator itPos = Q.begin();
 for(; itPos < Q.end(); itPos++) 
   {
     G4cout<<"----> Charge state "     <<(*itPos)->GetCharge()<<G4endl;
     G4cout<<"   => UnReactedFraction "<<(*itPos)->GetUnReactedFraction()<<G4endl;
     G4cout<<"   =>   ReactedFraction "<<(*itPos)->GetReactedFraction()<<G4endl;
     if((*itPos)->GetUseSetKEu())
       {
	 G4cout<<"   =>              KE/A "<<(*itPos)->GetSetKEu()/MeV<<" MeV"<<G4endl; 
	 //	 G4cout<<"   =>                KE "<<(*itPos)->GetSetKE()/MeV<<" MeV"<<G4endl; LR (not initialized if KEu is used)
       }
     else
       G4cout<<"   =>                KE "<<(*itPos)->GetSetKE()/MeV<<" MeV"<<G4endl;
   }
 CalcQUR();
 CalcQR();
}


//---------------------------------------------------------
void Outgoing_Beam::setDA(G4int da)
{

  DA=da;
 G4cout<<"----> Delta A for the outgoing beam set to  "<<DA<< G4endl; 
}
//---------------------------------------------------------
void Outgoing_Beam::setDZ(G4int dz)
{

  DZ=dz;
  G4cout<<"----> Delta Z for the outgoing beam set to  "<<DZ<< G4endl;
  
}
//---------------------------------------------------------
void Outgoing_Beam::setEx(G4double ex)
{

  Ex=ex;
  
G4cout<<"----> Excitation energy of the outgoing beam set to "<<
 G4BestUnit(Ex,"Energy")<<G4endl;

}
//---------------------------------------------------------
void Outgoing_Beam::setTarEx(G4double ex)
{

  TarEx=ex;
  
G4cout<<"----> Target excitation energy set to "<< G4BestUnit(TarEx,"Energy")<<G4endl;

}
//---------------------------------------------------------
void Outgoing_Beam::setTFrac(G4double ex)
{
  TFrac=ex;
  if(TFrac>1) TFrac=1.;
  if(TFrac<0) TFrac=0.;

  
G4cout<<"----> Fraction of target excitations set to "<<TFrac<<G4endl;

}
//-----------------------------------------------------------------
void Outgoing_Beam::settau(G4double t)
{

  tau=t;

G4cout<<"----> Lifetime of the excited state for the outgoing beam set to "<<
 G4BestUnit(tau,"Time")<<G4endl; 
}
//---------------------------------------------------------
void Outgoing_Beam::setbeta(G4double b)
{
// this function is defunct. was used for doppler reconstruction 
// when it was included in the simulation. Now it is handled separately
// TB
  betaDopp=b;

//G4cout<<"----> Beta for the outgoing beam set to "<<betaDopp<<G4endl; 
}
//---------------------------------------------------------
void Outgoing_Beam::setDoppZ(G4double z)
{
// TB same as setBeta
  posDopp.setZ(z);

//G4cout<<"----> Z position for the outgoing beam Doppler corrections set to "<<posDopp.getZ()<<G4endl; 
}
//---------------------------------------------------------
void Outgoing_Beam::setDoppY(G4double y)
{
// TB same as setbeta()
  posDopp.setY(y);

//G4cout<<"----> Y position for the outgoing beam Doppler corrections set to "<<posDopp.getY()<<G4endl; 
}

//---------------------------------------------------------
void Outgoing_Beam::setDoppX(G4double x)
{
// sameas setBeta()
  posDopp.setX(x);

//G4cout<<"----> X position for the outgoing beam Doppler corrections set to "<<posDopp.getX()<<G4endl; 
}
//---------------------------------------------------------
void Outgoing_Beam::SetUpChargeStates()
{

 vector<Charge_State*>::iterator itPos = Q.begin();
  // clear all elements from the array
  for(; itPos < Q.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  Q.clear();

  for(G4int i=0;i<NQ;i++) 
    Q.push_back(new Charge_State);


}
//---------------------------------------------------------
void Outgoing_Beam::SetQCharge(G4int q)
{

  if(SQ>=0&&SQ<NQ)
    {
      vector<Charge_State*>::iterator itPos = Q.begin();

      for(G4int i=0;i<SQ;i++) itPos++;
      (*itPos)->SetCharge(q);
    }
  else
    G4cout<<" Number of defined charge states ="<<NQ<<" is too small"<<G4endl;
}

//---------------------------------------------------------
void Outgoing_Beam::SetQUnReactedFraction(G4double f)
{

  if(SQ>=0&&SQ<NQ)
    {
      vector<Charge_State*>::iterator itPos = Q.begin();

      for(G4int i=0;i<SQ;i++) itPos++;
      (*itPos)->SetUnReactedFraction(f);
      CalcQUR();
    }
  else
    G4cout<<" Number of defined charge states ="<<NQ<<" is too small"<<G4endl;
}
//---------------------------------------------------------
void Outgoing_Beam::SetQReactedFraction(G4double f)
{

  if(SQ>=0&&SQ<NQ)
    {
      vector<Charge_State*>::iterator itPos = Q.begin();

      for(G4int i=0;i<SQ;i++) itPos++;
      (*itPos)->SetReactedFraction(f);
      CalcQR();
    }
  else
    G4cout<<" Number of defined charge states ="<<NQ<<" is too small"<<G4endl;
}
//---------------------------------------------------------
void Outgoing_Beam::SetQKEu(G4double e)
{
  G4int A;

  A=beamIn->getA()+DA;
  if(SQ>=0&&SQ<NQ)
    {
      vector<Charge_State*>::iterator itPos = Q.begin();

      for(G4int i=0;i<SQ;i++) itPos++;
      (*itPos)->SetKEu(e,A);
      CalcQR();
      CalcQUR();
    }
  else
    G4cout<<" Number of defined charge states ="<<NQ<<" is too small"<<G4endl;
}
//---------------------------------------------------------
void Outgoing_Beam::SetQKE(G4double e)
{

  if(SQ>=0&&SQ<NQ)
    {
      vector<Charge_State*>::iterator itPos = Q.begin();

      for(G4int i=0;i<SQ;i++) itPos++;
      (*itPos)->SetKE(e);
      CalcQR();
      CalcQUR();
    }
  else
    G4cout<<" Number of defined charge states ="<<NQ<<" is too small"<<G4endl;
}
//---------------------------------------------------------
void Outgoing_Beam::CalcQR()
{
  G4double sum;
  vector<Charge_State*>::iterator itPos = Q.begin();

  sum=0.;
  for(; itPos < Q.end(); itPos++)
    sum+=(*itPos)->GetReactedFraction();
  //  G4cout<<" Sum = "<<sum<<G4endl;
  itPos = Q.begin();
  G4int ind=0,i,max,n=0;
  for(; itPos < Q.end(); itPos++)
    { 
  
      max=(G4int)(1000.*(*itPos)->GetReactedFraction()/sum);
      //      G4cout<<" Sum = "<<sum<<" Max = "<<max<<" n "<<n<<G4endl;
      for(i=0;i<max;i++) 
	{
	  QR[ind]=(*itPos)->GetSetKE();
	  QRI[ind]=n;
	  ind++;
	}
      n++;
    }

}
//---------------------------------------------------------
void Outgoing_Beam::CalcQUR()
{
  G4double sum;
  vector<Charge_State*>::iterator itPos = Q.begin();

  sum=0.;
  for(; itPos < Q.end(); itPos++)
    sum+=(*itPos)->GetUnReactedFraction();

  itPos = Q.begin();
  G4int ind=0,i,max,n=0;
  for(; itPos < Q.end(); itPos++)
    { 
      max=(G4int)(1000.*(*itPos)->GetUnReactedFraction()/sum);
      //      G4cout<<" Sum = "<<sum<<" Max = "<<max<<" n "<<n<<" ind "<<ind<<G4endl;
      for(i=0;i<max;i++) 
	{
	  QUR[ind]=(*itPos)->GetSetKE();
	  QURI[ind]=n;
	  ind++;
	}
      n++;
    }
}
//---------------------------------------------------------
G4double Outgoing_Beam::GetURsetKE()
{
  G4int i;

  i=(G4int)(1000.*G4UniformRand());
  Index=QURI[i];
  return QUR[i];
}
//---------------------------------------------------------
G4double Outgoing_Beam::GetRsetKE()
{
  G4int i;

  i=(G4int)(1000.*G4UniformRand());
  Index=QRI[i];
  return QR[i];
}
//----------------------------------------------------------
void Outgoing_Beam::SetCoeff(int index,double a) {
  theAngularDistribution.SetCoeff(index,a);
  return;
}
//----------------------------------------------------------
void Outgoing_Beam::SetTargetCoeff(int index, double a)
{
  theTargetAngularDistribution.SetCoeff(index, a);
  return;
}
//----------------------------------------------------------
void Outgoing_Beam::openLvlSchemeFile()
{
  if(!lvlSchemeFile.is_open()) lvlSchemeFile.open(lvlSchemeFileName.c_str());
  
  if (lvlSchemeFile == NULL) G4cout << "lvlSchemeFile ERROR" << G4endl;
}
//----------------------------------------------------------
void Outgoing_Beam::closeLvlSchemeFile()
{
  lvlSchemeFile.close();

  return;
}
