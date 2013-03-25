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
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  particleTable = G4ParticleTable::GetParticleTable();
  BeamOut->SetReactionFlag(-1);

  // G4cout<<" +++++ In Generate Primaries "<<G4endl;
  if(source)
    {
      //      myDetector->GetSeGA()->DopplerOff();  //LR
      // G4cout<<" +++++ Shooting gammas "<<G4endl;
        particleGun->SetParticleDefinition(particleTable->FindParticle("gamma"));
     
        particleGun->SetParticleMomentumDirection(G4RandomDirection());
        particleGun->SetParticlePosition(sourcePosition);
	particleGun->SetParticleEnergy(GetSourceEnergy());

	//particleGun->SetParticlePosition(0);
    }
  if(inbeam)
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
    
	  //Reactions on the target
	  //TT=2.*myDetector->GetTarget()->GetZHalfLength();
	  TT= myDetector->GetTargetThickness();
	  TC=myDetector->GetTargetPlacement()->GetTranslation().getZ();
	  depth=TC+TT*(G4UniformRand()-0.5);
// 	      G4cout<< "- Target Thickness is  "<<TT/mm<<" mm"<<G4endl;
// 	      G4cout<< "- Target Center is at  "<<TC/mm<<" mm"<<G4endl;
// 	      G4cout<< "- Reaction depth   at  "<<depth/mm<<" mm"<<G4endl;
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
  } else if (name == "photopeaks") {
    SetSourcePhotopeaks();
  } else if (name == "eu152_peaks") {
    SetSourceEu152Peaks();
  } else if (name == "co56_peaks") {
    SetSourceCo56Peaks();
  } else if (name == "au") {
    SetSourceAu();
  } else if (name == "simple") {
    SetSourceSimple();
  }
  
}
//-------------------------------------------------------------
void PrimaryGeneratorAction::SetSourceEu152()
{
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
void PrimaryGeneratorAction::SetSourceCo60()
{
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

  rand=G4UniformRand()*sourceBranchingSum;

  vector<SourceData*>::iterator itPos = TheSource.begin();

  for(; itPos < TheSource.end(); itPos++)
    if(rand<(*itPos)->b) return (*itPos)->e;

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
