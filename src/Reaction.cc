#include "Reaction.hh"
#include "Reaction_Messenger.hh"

Reaction::Reaction(Outgoing_Beam* BO, const G4String& aName)
  : G4VProcess(aName), BeamOut(BO)
{
  if (verboseLevel>1) {
  G4cout <<GetProcessName() << " is created "<< G4endl;
     }
  BeamOut=BO;
  theProcessType=(G4ProcessType)6;       // Decay
  theProcessSubType=(G4ProcessType)231;  //DecayExt

  decayed_at_rest = false;
  target_reaction = false;
  ground_state    = false;

  //A cheap way to get a messenger
  theMessenger = new Reaction_Messenger(this);
}

Reaction::~Reaction() 
{
  delete theMessenger;
}

// In-flight reaction
// G4VParticleChange* Reaction::PostStepDoIt(
// 			     const G4Track& aTrack,
// 			     const G4Step& aStep
// 			    )
G4VParticleChange* Reaction::PostStepDoIt(
					  const G4Track& aTrack,
					  const G4Step&           // (unused parameter)
					  )
{

  // G4cout << "I'm in PostStepDoIt." << G4endl;
  // G4cout << "  " 
  // 	 << aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName()
  // 	 << G4endl;
  // G4cout << "  reaction_here: " 
  // 	 << reaction_here
  // 	 << G4endl;
  // G4cout << "  ground_state: " 
  // 	 << ground_state
  // 	 << G4endl;

  aParticleChange.Initialize(aTrack);

  if(reaction_here)
    {
      reaction_here=false;

      // Kill the track if we've already reacted and wandered back
      if(BeamOut->GetReactionFlag() == 1){

	aParticleChange.ProposeTrackStatus(fStopAndKill);

	// G4cout << "************************* PostStepDoIt: terminating track in "
	//        << aStep.GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
	//        << " at reaction depth"
	//        << G4endl;
       
      } 
      // React!
      else {

	//	G4cout << "*** PostStepDoIt: I'm reacting." << G4endl;
	
	BeamOut->ScanInitialConditions(aTrack);

	aParticleChange.ProposeTrackStatus(fStopAndKill);

	G4DynamicParticle *theProduct = NULL;
	if(BeamOut->AboveThreshold()){
	  aParticleChange.SetNumberOfSecondaries(1);
	  theProduct = BeamOut->ReactionProduct();
	  aParticleChange.AddSecondary(theProduct,
				       BeamOut->ReactionPosition(),
				       true);
	}

	BeamOut->SetReactionFlag(1);

#ifdef POL
	if(theProduct){
	  //This is a bit heavy, but lets me get the excitation energy
	  G4int Z = theProduct->GetParticleDefinition()->GetAtomicNumber();
	  G4int A = theProduct->GetParticleDefinition()->GetAtomicMass();
	  G4Fragment frag(A,Z,theProduct->Get4Momentum());
	  G4double Ex = frag.GetExcitationEnergy();
	  //Based on my reading of the source code, this should never return an
	  //invalid pointer...
	  const G4LevelManager *lman = G4NuclearLevelData::GetInstance()->GetLevelManager(Z,A);
	  size_t index = 0;
	  index = lman->NearestLevelIndex(Ex,index);
	  G4int twoJ = lman->SpinTwo(index);
	  if(substates.size() != (size_t)twoJ+1){
	    G4cerr << "Error: Reaction population parameter count != 2J+1."
		   << "\n  Use /reaction/population macro-file commands."
		   << "\n  (required when UCGretina is compiled with POL=1)"
		   << G4endl;
	    exit(EXIT_FAILURE);
	  }
	  std::vector<std::vector<G4complex>> polV;
	  polV.resize(twoJ+1);
	  G4double cosTheta = aTrack.GetMomentumDirection().cosTheta();
	  for(size_t k=0;k<polV.size();k++){
	    G4double PkKappa;
	    G4double Pk0_z = 0.;//First calculate along z-axis
	    for(int twoM=-twoJ;twoM<=twoJ;twoM+=2)
	      Pk0_z+=G4Clebsch::ClebschGordanCoeff(twoJ,twoM,twoJ,-twoM,2*k)/G4Clebsch::ClebschGordanCoeff(twoJ,twoM,twoJ,-twoM,0)/sqrt(2*k+1)*substates[twoM];
	    //Now rotate into particle frame
	    for(size_t kappa=0;kappa<k+1;kappa++){
	      PkKappa = Pk0_z*G4Clebsch::WignerLittleD(2*k,2*kappa,0,cosTheta);
	      polV[k].push_back(G4complex(PkKappa,0.0));
	    }
	  }
	  auto nucPstore = G4NuclearPolarizationStore::GetInstance();
	  auto pol = nucPstore->FindOrBuild(frag.GetZ_asInt(),frag.GetA_asInt(),frag.GetExcitationEnergy());
	  pol->SetPolarization(polV);
	}
#endif

      }

    }

  // Stop and kill the reaction product in its ground state.
  if(ground_state){
    ground_state = false;
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    // G4cout << "************************* PostStepDoIt: terminating track in "
    // 	   << aStep.GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
    // 	   << G4endl;

  }

  return &aParticleChange;
}

void Reaction::SetPopulation(G4int twoM, G4double frac){
  substates[twoM]=frac;
}

// Trigger the in-flight reaction at the depth in the target determined
// by the PrimaryGeneratorAction.
G4double Reaction::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
                             G4double,
                             G4ForceCondition* condition
                            )
{

  // G4cout << "I'm in PostStepGPIL." << G4endl;
  // G4cout << "  " 
  // 	 << aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName()
  // 	 << G4endl;
  // G4cout << "  momentum: " 
  // 	 << aTrack.GetDynamicParticle()->GetMomentum()
  // 	 << G4endl;
  // G4cout << "  target_reaction: " 
  // 	 << target_reaction
  // 	 << G4endl;

  reaction_here=false;
  *condition=NotForced;

  G4String name=aTrack.GetVolume()->GetLogicalVolume()->GetName();
  G4UserLimits* pUserLimits
    = aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();

  if(name=="Target_log"){
      
    // Target excitations:
    // Stop and kill the decay product once it reaches its ground state.
    if( target_reaction &&
	!G4StrUtil::contains(aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName(),'[') && BeamOut->GetReactionFlag()==1){
      ground_state = true;
      target_reaction = false;  //Reset for next decay
      return 0;
    }

    G4double ZReaction=pUserLimits->GetUserMinRange(aTrack);
    G4double ZCurrent=aTrack.GetPosition().getZ();
    G4double Z=ZReaction-ZCurrent;
    if(Z<0){
      // G4cout<<" Past the reaction point"<<G4endl;
      // G4cout<<" Volume "<<name<<G4endl;
      // G4cout<<" Z[mm]: reaction "<<ZReaction/mm<<" current "<<ZCurrent/mm<<" DZ "<<Z/mm<<G4endl;
      return DBL_MAX;
    } else if(Z>eps) {
      G4ThreeVector dir=aTrack.GetDynamicParticle()->GetMomentumDirection();
      
      dir*=(ZReaction-ZCurrent);
      // G4cout<<" Before the reaction point"<<G4endl;
      // G4cout<<" Volume "<<name<<G4endl;
      // G4cout<<" Z[mm]: reaction "<<ZReaction/mm<<" current "<<ZCurrent/mm<<" DZ "<<Z/mm<<G4endl;
      return dir.mag();
    } else if(Z<eps) {
      // G4cout<<" At the reaction point"<<G4endl;
      // G4cout<<" Volume "<<name<<G4endl;
      // G4cout<<" Z[mm]: reaction "<<ZReaction/mm<<" current "<<ZCurrent/mm<<" DZ "<<Z/mm<<G4endl;

      reaction_here = true;
      if( BeamOut->TargetExcitation() ) 
	target_reaction = true;
      return 0.;
    }
      
  }
    
  return DBL_MAX;
}

// If the reaction product comes to rest after emitting its gamma(s), kill it.
// G4VParticleChange* Reaction::AtRestDoIt(
// 			     const G4Track& aTrack,
// 			     const G4Step& aStep
// 			    )
G4VParticleChange* Reaction::AtRestDoIt(
			     const G4Track& aTrack,
			     const G4Step&          // unused parameter
			    )
{

  // G4cout << "I'm in AtRestDoIt." << G4endl;
  // G4cout << "  " 
  // 	 << aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName()
  // 	 << G4endl;

  aParticleChange.Initialize(aTrack);
    
  if( !G4StrUtil::contains(aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName(), '[') ) { 

    // G4cout << "************************* AtRestDoIt: terminating track in "
    // 	   << aStep.GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetName()
    // 	   << G4endl;

    aParticleChange.ProposeTrackStatus(fStopAndKill);

  }

  return &aParticleChange;

}
