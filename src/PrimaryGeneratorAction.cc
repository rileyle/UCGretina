#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction *detector, Incoming_Beam* BI, Outgoing_Beam* BO): BeamIn(BI), BeamOut(BO), myDetector(detector)
{
  SetInBeam();
  sourcePosition.setX(0);
  sourcePosition.setY(0);
#ifdef SCANNING
  sourcePosition.setY(200.*mm);
#endif
  sourcePosition.setZ(0);
  sourceRadius = 0;
  sourceDX = -1;
  sourceDY = -1;
  sourceSigmaX = -1;
  sourceSigmaY = -1;
  SetSourceEu152();
  sourceType = "eu152";
  n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  fracOn=false;
  frac=0;
  sourceWhiteLoE = -1;
  sourceWhiteHiE = -1;
  sourceMultiplicity = 1;
  isCollimated         = false;
  collimationDirection = G4ThreeVector(0., 0., 1.);
  collimationAngle     = 0.;
  thetaFileName = "";
  theta_max=180.*deg;
  theta_min=0.;
  theta_bin=0.;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  EventAction* eventAction
    = (EventAction*)G4RunManager::GetRunManager()->GetUserEventAction();
  if(eventAction->CacheIn())
    SetCache();
  
  particleTable = G4ParticleTable::GetParticleTable();
  ionTable = G4IonTable::GetIonTable();
  BeamOut->SetReactionFlag(-1);
  
  if(source)
    {

      if(sourceType == "muon")
	{
	  if(G4UniformRand()<0.5)
	    particleGun->SetParticleDefinition(particleTable->FindParticle("mu+"));
	  else
	    particleGun->SetParticleDefinition(particleTable->FindParticle("mu-"));
	}
      else if(sourceType == "neutron")
	{
	  particleGun->SetParticleDefinition(particleTable->FindParticle("neutron"));
	}
      else
	{
	  particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
	}
      
      if(background)
	{
	  if(sourceType == "muon") // Cosmic-ray background
	    {
	      // Cosmic rays with a uniform distribution on a
	      // circle with diameter covering the crystals. 
	      // Maybe.
	      G4double x = 28.*cm*(1.0 - 2*G4UniformRand());
	      G4double z = sqrt(28.*cm*28.*cm - x*x)*(1.0 - 2*G4UniformRand());
	      particleGun->SetParticlePosition(G4ThreeVector(x, myDetector->GetBGSphereRmin(), z));
	      // Vertical for now. Eventually implement cos^2 distribution
	      particleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1,0));
	    }
	  else // Room background
	    {
	      G4ThreeVector v1, v2;
	      // Inside the background sphere walls.
	      G4double bgRmin = myDetector->GetBGSphereRmin();
	      G4double bgRmax = myDetector->GetBGSphereRmax();
	      v1 = G4RandomDirection()*(bgRmin + 20.0*cm +
					G4UniformRand()*(bgRmax -
							 bgRmin - 20.0*cm));  
	      // Inside the outer surface of the Gretina mounting sphere
	      v2 = G4RandomDirection()*G4UniformRand()*1.276/2.*m;  
	      particleGun->SetParticlePosition(v1);
	      particleGun->SetParticleMomentumDirection(v2-v1);
	    }
	}
      else
	{
	  // Find the emission point.
	  G4ThreeVector sourceOffset = G4ThreeVector(0.,0.,0.);
	  G4double xOffset = 0;
	  G4double yOffset = 0;
	  if(sourceDX > 0 || sourceDY > 0 ||
	     sourceSigmaX > 0 || sourceSigmaY > 0 ){
	    if(sourceDX > 0)
	      xOffset = (0.5-G4UniformRand())*sourceDX;
	    else if(sourceSigmaX > 0)
	      xOffset = CLHEP::RandGauss::shoot(0, sourceSigmaX);
	    if(sourceDY > 0)
	      yOffset = (0.5-G4UniformRand())*sourceDY;
	    else if(sourceSigmaY > 0)
	      yOffset = CLHEP::RandGauss::shoot(0, sourceSigmaY);
	  } else if(sourceRadius > 0){
	    G4double rOffset, phiOffset;
	    phiOffset = G4UniformRand()*8.*atan(1.);
	    rOffset   = G4UniformRand()+G4UniformRand();
	    if(rOffset >= 1) rOffset =- (rOffset - 2.);
	    xOffset = rOffset*cos(phiOffset)*sourceRadius;
	    yOffset = rOffset*sin(phiOffset)*sourceRadius;
	  }
	  sourceOffset.setX(xOffset);
	  sourceOffset.setY(yOffset);

	  if(isCollimated){
	    G4double cosThetaMax = cos(collimationAngle);
	    G4double dcosTheta = 1 - cos(collimationAngle);
	    G4double cosTheta =  cosThetaMax + G4UniformRand()*dcosTheta;
	    G4double sinTheta2 = 1. - cosTheta*cosTheta;
	    if( sinTheta2 < 0.)  sinTheta2 = 0.;
	    G4double sinTheta  = std::sqrt(sinTheta2); 
	    G4double phi       = twopi*G4UniformRand();
	    G4ThreeVector sourceDir = G4ThreeVector(sinTheta*std::cos(phi),
						    sinTheta*std::sin(phi), 
						    cosTheta).unit();
	    sourceDir.rotateY(collimationDirection.theta());
	    sourceDir.rotateZ(collimationDirection.phi());
	    particleGun->SetParticleMomentumDirection(sourceDir);
	    if(sourceRadius > 0){
	      sourceOffset.rotateY(collimationDirection.theta());
	      sourceOffset.rotateZ(collimationDirection.phi());
	    }
	  } else {
	    if(thetaFileName == ""){
	      particleGun->SetParticleMomentumDirection(G4RandomDirection());
	    } else {
	      CLHEP::RandGeneral randTheta(thetaDist, Ntheta);
	      G4double theta
		= theta_min + randTheta.fire()*(theta_max - theta_min);
	      G4double phi = G4UniformRand()*twopi;
	      
	      G4ThreeVector sourceDir
		= G4ThreeVector(std::sin(theta)*std::cos(phi),
				std::sin(theta)*std::sin(phi), 
				std::cos(theta)).unit();

	      particleGun->SetParticleMomentumDirection(sourceDir);
	    }
	  }
	  particleGun->SetParticlePosition(sourcePosition+sourceOffset);
	}

      particleGun->SetParticleEnergy(GetSourceEnergy());

      //  G4cout<<" +++++ Generating an event "<<G4endl;
      particleGun->GeneratePrimaryVertex(anEvent);
      PrimaryVertexInformation* primaryVertexInfo = new PrimaryVertexInformation;
      anEvent->GetPrimaryVertex()->SetUserInformation(primaryVertexInfo);

    }
  else if(inbeam)
    {

      ion =  ionTable->GetIon(BeamIn->getZ(),BeamIn->getA(),BeamIn->getEx());
      particleGun->SetParticleDefinition(ion);
      
      direction=BeamIn->getDirection();
      particleGun->SetParticleMomentumDirection(direction);

      position=BeamIn->getPosition(direction);

      // Include target offset (if any).
      G4ThreeVector offset;
      offset.setX(myDetector->GetTargetPos()->x());
      offset.setY(myDetector->GetTargetPos()->y());
      offset.setZ(myDetector->GetTargetPos()->z());
      position += offset;

      particleGun->SetParticlePosition(position);
      
      KE=BeamIn->getKE(ion);
      particleGun->SetParticleEnergy(KE);

      G4double TT;
      G4double TC;
      G4double depth;
    
      //Reactions in the target
#ifdef LHTARGET
      TC = myDetector->GetTargetPos()->getZ(); // Z position of the target center

      // Intersection of the beam trajectory with the central plane of the target
      G4ThreeVector mid = position
        + direction*(TC - position.getZ())/direction.getZ(); 

      // The DistanceToOut() method works in the reference frame of
      // the logical volume, not the physical one. The origin of
      // this frame for the 200 mg cell is at the center of the upstream
      // surface, so the mid vector is shifted accordingly.
      if(myDetector->GetTargetThickness() >= 29.0*mm)
	mid.setZ(myDetector->GetTargetThickness()/2.0);
      else
	mid.setZ(0);
            
      // Distance along current trajectory in the target volume
      TT = myDetector->GetTarget()->DistanceToOut(mid, direction)
	+ myDetector->GetTarget()->DistanceToOut(mid, -direction);
#else
      TT = myDetector->GetTarget()->DistanceToIn(position, direction);
      TT = myDetector->GetTarget()->DistanceToOut(position+TT*direction, 
						  direction);

      TC=myDetector->GetTargetPos()->getZ();
#endif
      TT *= direction.getZ();
      
      // Z position of the reaction point
      depth=TC+TT*(G4UniformRand()-0.5);
      myDetector->setTargetReactionDepth(depth);

      // G4cout << "**** center = " << TC
      // 	     << ", direction = " << direction
      // 	     << ", thickness along trajectory = " << TT
      // 	     << ", depth = " << depth
      // 	     << G4endl;

      //  G4cout<<" +++++ Generating an event "<<G4endl;
      particleGun->GeneratePrimaryVertex(anEvent);
      PrimaryVertexInformation* primaryVertexInfo
	= new PrimaryVertexInformation;
      anEvent->GetPrimaryVertex()->SetUserInformation(primaryVertexInfo);

    }
  else if(cache)
    {
      G4int Npoints;
      G4double x[1000];
      G4double y[1000];
      G4double z[1000];
      G4double bx[1000];
      G4double by[1000];
      G4double bz[1000];
      G4double t[1000];
      G4double ata, bta, dta, yta;

      // Use the half life ( /Experiment/Cache/HalfLife ) to find the
      // emission position.
      G4bool emissionSet = false;
      G4ThreeVector emissionPos;
      G4ThreeVector emissionBeta;

#ifdef CACHETEXT
      char currentLine[1000];
      std::ifstream& cacheFile = eventAction->getCacheInputFile();
      cacheFile >> Npoints;
      if(cacheFile.eof()) {
	G4cerr << "Error reading trajectory header from cache file."
	       << G4endl;
	exit(EXIT_FAILURE);
      }
      //	G4cout << Npoints << G4endl;
#else
      std::FILE* cacheFile = eventAction->getCacheInputFile();
      G4int size = fread(&Npoints, sizeof(G4int), 1, cacheFile);
      if(size != 1) {
	G4cerr << "Error reading trajectory header from cache file."
	       << G4endl;
	exit(EXIT_FAILURE);
      }
#endif
      G4double halflife = eventAction->GetCacheHalfLife(); // ps
      G4double lifetime = halflife/log(2)/ps;
      G4double emissionTime = CLHEP::RandExponential::shoot(lifetime);
      for(int i = 0; i < Npoints; i++) {
#ifdef CACHETEXT
	if(cacheFile.eof()) {
	  G4cerr << "Error reading trajectory point from cache file."
		 << G4endl;
	  exit(EXIT_FAILURE);
	}
	cacheFile >> x[i] >> y[i] >> z[i] >> bx[i] >> by[i] >> bz[i] >> t[i];
	// G4cout << x[i] << std::setw(12)
	// 	 << y[i] << std::setw(12)
	// 	 << z[i] << std::setw(12)
	// 	 << bx[i] << std::setw(12)
	// 	 << by[i] << std::setw(12)
	// 	 << bz[i] << std::setw(12)
	// 	 << t[i] << std::setw(12) << G4endl;
	cacheFile.getline(currentLine,1000);
#else
	CTP inputPosData;
	size = fread(&inputPosData, sizeof(inputPosData), 1, cacheFile);
	if(size != 1) {
	  G4cerr << "Error reading trajectory point from cache file."
		 << G4endl;
	  exit(EXIT_FAILURE);
	}
	x[i] = inputPosData.x;
	y[i] = inputPosData.y;
	z[i] = inputPosData.z;
	bx[i] = inputPosData.bx;
	by[i] = inputPosData.by;
	bz[i] = inputPosData.bz;
	t[i] = inputPosData.t;
#endif
	if(!emissionSet){
	  if(lifetime > 0){
	    // Inside the target, we interpolate.
	    if(t[i] > emissionTime){
	      G4ThreeVector x1 = G4ThreeVector(x[i-1], y[i-1], z[i-1]);
	      G4ThreeVector x2 = G4ThreeVector(x[i],   y[i],   z[i]);
	      emissionPos = (x1*(emissionTime-t[i-1])
			     + x2*(t[i]-emissionTime))/(t[i] - t[i-1]);
	      G4ThreeVector beta1 = G4ThreeVector(bx[i-1], by[i-1], bz[i-1]);
	      G4ThreeVector beta2 = G4ThreeVector(bx[i],   by[i],   bz[i]);
	      emissionBeta = (beta1*(emissionTime-t[i-1])
			      + beta2*(t[i]-emissionTime))/(t[i] - t[i-1]);
	      emissionSet = true;
	    }
	  } else {
	    emissionPos = G4ThreeVector(x[0],   y[0],  z[0]);
	    emissionBeta = G4ThreeVector(bx[0], by[0], bz[0]);
	    emissionSet = true;
	  }
	}
      }
      // Downstream of the target, we extrapolate.
      if(!emissionSet){
	emissionBeta = G4ThreeVector(bx[Npoints-1], by[Npoints-1], bz[Npoints-1]);
	G4double c = 0.299792458; // mm/ps
	emissionPos = G4ThreeVector(x[Npoints-1], y[Npoints-1], z[Npoints-1])
	  + emissionBeta*c*(emissionTime-t[Npoints-1]);
	emissionSet = true;
      }
      
      // Read S800 data.
#ifdef CACHETEXT
      if(cacheFile.eof()) {
	G4cerr << "Error reading S800 data from cache file."
	       << G4endl;
	exit(EXIT_FAILURE);
      }
      cacheFile >> ata >> bta >> dta >> yta;
      cacheFile.getline(currentLine,1000);
      //	G4cout << ata << std::setw(12) << bta << std::setw(12) << dta << std::setw(12) << yta << std::setw(12) << G4endl;
#else
      CS inputS8Data;
      size =  fread(&inputS8Data, sizeof(inputS8Data), 1, cacheFile);
      if(size != 1) {
	G4cerr << "Error reading S800 data from cache file."
	       << G4endl;
	exit(EXIT_FAILURE);
      }
      ata = inputS8Data.ata;
      bta = inputS8Data.bta;
      dta = inputS8Data.dta;
      yta = inputS8Data.yta;
#endif
      //G4FourVector r = (ata, bta, dta, yta);
      // Projectile-frame gamma-ray energy
      G4double e = eventAction->GetCacheGammaEnergy();

      // Projectile-frame gamma-ray momentum
      G4ThreeVector p;
      AngularDistribution* angDist = eventAction->GetCacheAngDist();
      if(angDist == NULL){ // Isotropic
	p = G4RandomDirection();
      } else {             // User-specified a0, a2, a4
	G4double theta = angDist->GetRandomAngle();
	G4double cosTheta = std::cos(theta);
	G4double sinTheta = std::sin(theta);
	G4double phi = twopi * G4UniformRand();
	p.set(sinTheta*std::cos(phi),
	      sinTheta*std::sin(phi),
	      cosTheta);
      }
      p.setMag(e);
      
      // Projectile-frame gamma-ray 4 momentum
      G4LorentzVector pProj = G4LorentzVector(p, e);

      // Boost the gamma-ray 4 momentum to the lab frame.
      G4LorentzVector pLab = pProj.boost(emissionBeta);
      
      // Lab-frame gamma-ray emission direction
      G4ThreeVector gDir = G4ThreeVector(0, 0, 1); 
      gDir.setTheta(pLab.theta());
      gDir.setPhi(pLab.phi());

      //      G4cout << "emission Pos = " << emissionPos
      //	     << ", emission Beta = " << emissionBeta
      //	     << ", emission Time = " << emissionTime
      //	     << ", gDir = " << gDir
      //	     << ", pLab.e = " << pLab.e()
      //	     << G4endl;
      
      // Set the particle gun parameters.
      particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
      particleGun->SetParticleMomentumDirection(gDir);
      particleGun->SetParticleEnergy(pLab.e());
      particleGun->SetParticlePosition(emissionPos);
      particleGun->GeneratePrimaryVertex(anEvent);

      PrimaryVertexInformation* primaryVertexInfo = new PrimaryVertexInformation;
      anEvent->GetPrimaryVertex()->SetUserInformation(primaryVertexInfo);
      primaryVertexInfo->AddBeta(emissionBeta.mag(), 0);
      primaryVertexInfo->SetATA(ata);
      primaryVertexInfo->SetBTA(bta);
      primaryVertexInfo->SetDTA(dta);
      primaryVertexInfo->SetYTA(yta);
  
    }

  if(sourceMultiplicity > 1 &&
     (sourceType == "white" || sourceType == "bgwhite")){
    for(G4int i = 1; i<sourceMultiplicity; i++){
      particleGun->SetParticleEnergy(GetSourceEnergy());
      particleGun->GeneratePrimaryVertex(anEvent);
    }
  }
  // G4cout<<" +++++ Out Generate Primaries "<<G4endl;
}
//---------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceOnTargetFace()
{
  sourcePosition.setX(0);
  sourcePosition.setY(0);
  //sourcePosition.setZ(-myDetector->GetTarget()->GetZHalfLength()+ //LR for LH
  sourcePosition.setZ(-myDetector->GetTargetThickness()/2.0+
	myDetector->GetTargetPlacement()->GetTranslation().getZ());
  particleGun->SetParticlePosition(sourcePosition);
}

//---------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceOnTargetBack()
{
  sourcePosition.setX(0);
  sourcePosition.setY(0);
  //sourcePosition.setZ(myDetector->GetTarget()->GetZHalfLength()+ //LR for LH
  sourcePosition.setZ(myDetector->GetTargetThickness()/2.0+
	myDetector->GetTargetPlacement()->GetTranslation().getZ());
  particleGun->SetParticlePosition(sourcePosition);
}

//---------------------------------------------------------------------
void PrimaryGeneratorAction::SourceReport()
{
  if(source)
    {
      G4cout<<"----> Source type is set to "<< sourceType << G4endl; //LR
      G4cout<<"----> Source position in X is set to "<<G4BestUnit(sourcePosition.getX(),"Length")<<G4endl;
      G4cout<<"----> Source position in Y is set to "<<G4BestUnit(sourcePosition.getY(),"Length")<<G4endl;
      G4cout<<"----> Source position in Z is set to "<<G4BestUnit(sourcePosition.getZ(),"Length")<<G4endl;

      if(sourceDX > 0 || sourceDY > 0 ||
	 sourceSigmaX > 0 || sourceSigmaY > 0 ){
	if(sourceDX > 0)
	  G4cout<<"----> Source width (nondispersive) is set to "<<G4BestUnit(sourceDX,"Length")<<G4endl;
	if(sourceDY > 0)
	  G4cout<<"----> Source height (dispersive) is set to "<<G4BestUnit(sourceDX,"Length")<<G4endl;
	if(sourceSigmaX > 0)
	  G4cout<<"----> Source Gaussian width (nondispersive) is set to "<<G4BestUnit(sourceSigmaX,"Length")<<G4endl;
	if(sourceSigmaY > 0)
	  G4cout<<"----> Source Gaussian height (dispersive) is set to "<<G4BestUnit(sourceSigmaY,"Length")<<G4endl;

      } else {
	G4cout<<"----> Source disk radius is set to "<<G4BestUnit(sourceRadius,"Length")<<G4endl;
      }
    }
  else
    G4cout<<"----> In-beam run selected for simulations"<<G4endl;
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceType(G4String name) //LR
{

  sourceType = name;

  if(name == "eu152") {
    SetSourceEu152();
  } else if (name == "cs137") {
    SetSourceCs137();
  } else if (name == "co56") {
    SetSourceCo56();
  } else if (name == "co60") {
    SetSourceCo60();
  } else if (name == "ra226") {
    SetSourceRa226();
  } else if (name == "am241") {
    SetSourceAm241();
  } else if (name == "photopeaks") {
    SetSourcePhotopeaks();
  } else if (name == "eu152_peaks") {
    SetSourceEu152Peaks();
  } else if (name == "co56_peaks") {
    SetSourceCo56Peaks();
  } else if (name == "ra226_peaks") {
    SetSourceRa226Peaks();
  } else if (name == "au") {
    SetSourceAu();
  } else if (name == "white") {
    SetSourceWhite();
  } else if (name == "background") {
    SetSourceBG();
  } else if (name == "bgwhite") {
    SetSourceBGWhite();
  } else if (name == "muon") {
    SetSourceMuon();
  } else if (name == "neutron") {
    SetSourceNeutron();
  } else if (name == "simple") {
    SetSourceSimple();
  }
  
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceEu152()
{
  sourceType = "eu152";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=121.782*keV;sourceBranchingSum+=13607.;                  //LR  13620.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 244.699*keV; sourceBranchingSum+= 3612.;                //LR   3590.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 295.939*keV; sourceBranchingSum+=  211.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 344.281*keV; sourceBranchingSum+=12743.;                //LR  12750.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 367.789*keV; sourceBranchingSum+=  405.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 411.126*keV; sourceBranchingSum+= 1073.;                 //LR   1071.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 443.965*keV; sourceBranchingSum+= 1499.;                 //LR   1480.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 488.661*keV; sourceBranchingSum+=  195.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 564.021*keV; sourceBranchingSum+=  236.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 586.294*keV; sourceBranchingSum+=  220.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 678.578*keV; sourceBranchingSum+=  221.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 688.678*keV; sourceBranchingSum+=  400.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 778.903*keV; sourceBranchingSum+= 6221.;                 //LR   6190.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 867.390*keV; sourceBranchingSum+= 2021.;                 //LR   1990.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 964.055*keV; sourceBranchingSum+= 7017.;                 //LR   6920.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1005.279*keV; sourceBranchingSum+=  310.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1085.842*keV; sourceBranchingSum+= 4859.;                 //LR
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1089.767*keV; sourceBranchingSum+=  830.;                 //LR  1089.7, 820.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1109.180*keV; sourceBranchingSum+=   88.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1112.087*keV; sourceBranchingSum+= 6494.;                 //LR 6590.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1212.970*keV; sourceBranchingSum+=  677.;                 //LR 670.
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1299.152*keV; sourceBranchingSum+=  780.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1408.022*keV; sourceBranchingSum+=10000.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceEu152Peaks()
{
  sourceType = "eu152_peaks";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=121.782*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 244.699*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 344.281*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 411.126*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 443.965*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 778.903*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 964.055*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1112.087*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1408.022*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceCs137()
{
  sourceType = "cs137";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=661.657*keV;sourceBranchingSum+=100.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceCo56()
{
  sourceType = "co56";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=846.764*keV;sourceBranchingSum+=99.933;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=977.373*keV; sourceBranchingSum+= 1.449;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1037.884*keV; sourceBranchingSum+= 14.13;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1175.099*keV; sourceBranchingSum+=2.239;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1238.287*keV; sourceBranchingSum+=66.07;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1360.206*keV; sourceBranchingSum+=4.256;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1771.350*keV; sourceBranchingSum+=15.49;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2015.179*keV; sourceBranchingSum+=3.029;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2034.759*keV; sourceBranchingSum+=7.71;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2598.460*keV; sourceBranchingSum+=16.96;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3009.596*keV; sourceBranchingSum+=1.16;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3201.954*keV; sourceBranchingSum+=3.13;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3253.417*keV; sourceBranchingSum+=7.62;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3272.998*keV; sourceBranchingSum+=1.78;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3451.154*keV; sourceBranchingSum+=0.93;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3548.27*keV; sourceBranchingSum+=0.178;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceCo56Peaks()
{
  sourceType = "co56_peaks";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=846.764*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1037.884*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1175.099*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1238.287*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1360.206*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1771.350*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2034.759*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2598.460*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3201.954*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3451.154*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3548.27*keV; sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceRa226()
{
  sourceType = "ra226";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=186.211*keV;sourceBranchingSum+=0.03533;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=241.997*keV;sourceBranchingSum+=0.0719;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=295.224*keV;sourceBranchingSum+=0.1828;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=351.932*keV;sourceBranchingSum+=0.3534;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=609.316*keV;sourceBranchingSum+=0.4516;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=665.453*keV;sourceBranchingSum+=0.01521;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=768.367*keV;sourceBranchingSum+=0.0485;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=806.185*keV;sourceBranchingSum+=0.01255;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=934.061*keV;sourceBranchingSum+=0.03074;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1120.287*keV;sourceBranchingSum+=0.1478;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1155.190*keV;sourceBranchingSum+=0.01624;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1238.110*keV;sourceBranchingSum+=0.05785;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1280.960*keV;sourceBranchingSum+=0.01425;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1377.669*keV;sourceBranchingSum+=0.03954;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1401.516*keV;sourceBranchingSum+=0.01324;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1407.993*keV;sourceBranchingSum+=0.02369;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1509.217*keV;sourceBranchingSum+=0.02108;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1620.740*keV;sourceBranchingSum+=0.0151;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1661.316*keV;sourceBranchingSum+=0.01037;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1729.640*keV;sourceBranchingSum+=0.02817;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1764.539*keV;sourceBranchingSum+=0.1517;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1847.420*keV;sourceBranchingSum+=0.0200;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2118.536*keV;sourceBranchingSum+=0.01148;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2204.071*keV;sourceBranchingSum+=0.0489;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2447.673*keV;sourceBranchingSum+=0.01536;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));

}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceAm241()
{
  sourceType = "am241";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=26.3446*keV;sourceBranchingSum+=0.0240;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=33.1963*keV;sourceBranchingSum+=0.00121;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=59.5409*keV;sourceBranchingSum+=0.3578;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));

}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceRa226Peaks()
{
  sourceType = "ra226_peaks";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=186.211*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=241.997*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=295.224*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=351.932*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=609.316*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=768.367*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1120.287*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1764.539*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2204.071*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2447.673*keV;sourceBranchingSum+=1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));

}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceCo60()
{
  sourceType = "co60";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=1173.238*keV;sourceBranchingSum+=99.857;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1332.502*keV; sourceBranchingSum+=99.983;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourcePhotopeaks()
{
  sourceType = "photopeaks";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e= 121.782*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 244.699*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 344.281*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 411.126*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 443.965*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 778.903*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 964.055*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1112.087*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1408.022*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e= 846.764*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1037.884*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1175.099*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1238.287*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1360.206*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1771.350*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2034.759*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2598.460*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3201.954*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3451.154*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=3548.270*keV; sourceBranchingSum+= 1.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceAu()
{
  sourceType = "au";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=547.5*keV;sourceBranchingSum+=100.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceSimple()
{
  sourceType = "simple";

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=1000.0*keV;sourceBranchingSum+=100.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceWhite()
{
  sourceType = "white";

  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();
}
//-------------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceBGWhite()
{
  sourceType = "bgwhite";

  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();
}
//-------------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceBG()
{
  background = true;

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
  // finally, clear all elements from the array
  TheSource.clear();

  // 232Th, 238U relative intensities from:
  // Gordon R. Gilmore,
  // Practical Gamma-Ray Spectrometry, 2nd Edition, Wiley
  // Published Online: 21 APR 2008
  
  e=46.54*keV;sourceBranchingSum+=4.25;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=63.28*keV;sourceBranchingSum+=4.8;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=74.82*keV;sourceBranchingSum+=27.7;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=77.11*keV;sourceBranchingSum+=46.2;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=92.58*keV;sourceBranchingSum+=5.58;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=185.72*keV;sourceBranchingSum+=57.2;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=238.63*keV;sourceBranchingSum+=43.6;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=295.22*keV;sourceBranchingSum+=18.5;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=351.93*keV;sourceBranchingSum+=35.6;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=583.19*keV;sourceBranchingSum+=30.6;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=609.31*keV;sourceBranchingSum+=45.49;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=911.20*keV;sourceBranchingSum+=25.8;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=968.97*keV;sourceBranchingSum+=15.8;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1460.82*keV;sourceBranchingSum+=100.; // NSCL S3 Vault June 2013
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=1764.49*keV;sourceBranchingSum+=15.28;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
  e=2614.51*keV;sourceBranchingSum+=35.38;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));
}
//-------------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceMuon()
{
  sourceType = "muon";
  background = true;

  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  // Approximate average ground level muon energy according to the PDG.
  e=4.0*GeV;sourceBranchingSum+=100.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));

}
//-------------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceNeutron()
{
  sourceType = "neutron";
  background = false;

#ifndef NEUTRONS
  G4cout << "Error: To use the neutron source type, compile with NEUTRONS=1."
	 << G4endl;
  exit(EXIT_FAILURE);
#endif
  
  G4double e;
  sourceBranchingSum=0.;

  // start from the beginning of the array
  vector<SourceData*>::iterator itPos = TheSource.begin();
  // clear all elements from the array
  for(; itPos < TheSource.end(); itPos++)
    delete *itPos;    // free the element from memory
   // finally, clear all elements from the array
  TheSource.clear();

  e=1000.0*keV;sourceBranchingSum+=100.;
  TheSource.push_back(new SourceData(e,sourceBranchingSum));

}
//-------------------------------------------------------------------------
G4double PrimaryGeneratorAction::GetSourceEnergy()
{
 
  G4double rand;
  if(sourceWhiteLoE >= 0 && sourceWhiteHiE > sourceWhiteLoE){

    return sourceWhiteLoE + G4UniformRand()*(sourceWhiteHiE - sourceWhiteLoE);

  } else {

    rand=G4UniformRand()*sourceBranchingSum;

    vector<SourceData*>::iterator itPos = TheSource.begin();

    for(; itPos < TheSource.end(); itPos++)
      if(rand<(*itPos)->b) return (*itPos)->e;

  }

  cout << "******** Oops!!!!" << endl;

  return -1*keV;
}
//-------------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceEnergy(G4double energy)
{
  if(sourceType == "simple" || sourceType == "neutron"){
     vector<SourceData*>::iterator itPos = TheSource.begin();
     (*itPos)->e = energy;
     G4cout << "Setting source energy to " << energy << " MeV" << G4endl;
  } else {
    G4cout << "Warning: /Experiment/Source/setEnergy has no effect unless the source type is set to \"simple\" or  \"neutron\"" << G4endl;
  }
}
//-------------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceThetaFile(G4String fileName)
{
 
  char line[1000];

  thetaFileName = fileName;
  std::ifstream thetaFile;
  thetaFile.open(thetaFileName);

  thetaFile >> theta_min >> theta_max >> theta_bin;
  thetaFile.getline(line,1000);  // Advance to next line.

  G4double x;
  Ntheta = 0;
  while(thetaFile >> x){
    thetaFile.getline(line,1000);  // Advance to next line.
    thetaDist[Ntheta] = x;
    Ntheta++;
  }

  thetaFile.close();
  
}
//-------------------------------------------------------------------------
