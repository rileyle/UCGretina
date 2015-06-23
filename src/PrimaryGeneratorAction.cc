#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction *detector, Incoming_Beam* BI, Outgoing_Beam* BO): BeamIn(BI), BeamOut(BO), myDetector(detector)
{
  SetInBeam();
  sourcePosition.setX(0);
  sourcePosition.setY(0);
  sourcePosition.setZ(0);
  SetSourceEu152();
  sourceType = "eu152";
  n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  fracOn=false;
  frac=0;
  sourceWhiteLoE = 100.*keV;
  sourceWhiteHiE = 10000.*keV;
  isCollimated         = false;
  collimationDirection = G4ThreeVector(0., 0., 1.);
  collimationAngle     = 0.;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  particleTable = G4ParticleTable::GetParticleTable();
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
	  } else {
	    particleGun->SetParticleMomentumDirection(G4RandomDirection());
	  }
	  particleGun->SetParticlePosition(sourcePosition);
	}

      particleGun->SetParticleEnergy(GetSourceEnergy());

    }
  else if(inbeam)
    {

      ion =  particleTable->GetIon(BeamIn->getZ(),BeamIn->getA(),BeamIn->getEx());
      particleGun->SetParticleDefinition(ion);
      
      position=BeamIn->getPosition();
      particleGun->SetParticlePosition(position);
      
      direction=BeamIn->getDirection();
      particleGun->SetParticleMomentumDirection(direction);
      
      KE=BeamIn->getKE(ion);
      particleGun->SetParticleEnergy(KE);

      if(fracOn)
	{
	  if(G4UniformRand()<frac) 
	    BeamOut->SetReactionOff();
	  else   
	    BeamOut->SetReactionOn();
	}

      if(BeamOut->ReactionOn())
      	{

	  G4double TT;
	  G4double TC;
	  G4double depth;
    
	  //Reactions in the target

	  // This works in general ...
	  TT = myDetector->GetTarget()->DistanceToIn(position, direction);
	  TT = myDetector->GetTarget()->DistanceToOut(position+TT*direction, 
						      direction);
          TT *= direction.getZ();

	  // ... but this may be faster, and approximately correct for 
	  // a flat target.
	  //	  TT= myDetector->GetTargetThickness();

	  TC=myDetector->GetTargetPos()->getZ();
	  depth=TC+TT*(G4UniformRand()-0.5);

 	  //    G4cout<< "- Target Thickness is  "<<TT/mm<<" mm"<<G4endl;
 	  //    G4cout<< "- Target Center is at  "<<TC/mm<<" mm"<<G4endl;
 	  //    G4cout<< "- Reaction depth   at  "<<depth/mm<<" mm"<<G4endl;
 	  //    G4cout<< "- Direction is  "<<direction<<G4endl;

	  myDetector->setTargetReactionDepth(depth);

	}

    }

  // G4cout<<" +++++ Generating an event "<<G4endl;
  particleGun->GeneratePrimaryVertex(anEvent);
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
G4double PrimaryGeneratorAction::GetSourceEnergy()
{
 
  G4double rand;
  if(sourceType != "white" && sourceType != "bgwhite"){

    rand=G4UniformRand()*sourceBranchingSum;

    vector<SourceData*>::iterator itPos = TheSource.begin();

    for(; itPos < TheSource.end(); itPos++)
      if(rand<(*itPos)->b) return (*itPos)->e;

  } else {
    
    return sourceWhiteLoE + G4UniformRand()*(sourceWhiteHiE - sourceWhiteLoE);

  }

  cout << "******** Oops!!!!" << endl;

  return -1*keV;
}
//-------------------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceEnergy(G4double energy)
{
  if(sourceType == "simple"){
     vector<SourceData*>::iterator itPos = TheSource.begin();
     (*itPos)->e = energy;
     G4cout << "Setting source energy to " << energy << " MeV" << G4endl;
  } else {
    G4cout << "Warning: /Experiment/Source/setEnergy has no effect unless the source type is set to \"simple\"" << G4endl;
  }
}
//-------------------------------------------------------------------------
