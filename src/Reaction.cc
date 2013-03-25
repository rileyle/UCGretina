#include "Reaction.hh"



Reaction::Reaction(Outgoing_Beam* BO, const G4String& aName)
  : G4VProcess(aName), BeamOut(BO)
{
  if (verboseLevel>1) {
  G4cout <<GetProcessName() << " is created "<< G4endl;
     }
  BeamOut=BO;
  theProcessType=(G4ProcessType)6;       // Decay
  theProcessSubType=(G4ProcessType)231;  //DecayExt
}

Reaction::~Reaction() 
{                                     
}                                     

G4VParticleChange* Reaction::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step& 
			    )
//
// Stop the current particle, if requested by G4UserLimits 
// 			    			    			    
{

  aParticleChange.Initialize(aTrack);

  if(reaction_here)
    {
      reaction_here=false;

      BeamOut->ScanInitialConditions(aTrack);

      aParticleChange.ProposeTrackStatus(fStopAndKill);

 

      if(G4UniformRand()<BeamOut->getTFrac())
	{
	  aParticleChange.SetNumberOfSecondaries(2);
	  aParticleChange.AddSecondary(BeamOut->ProjectileGS(),BeamOut->ReactionPosition(),true);
	  aParticleChange.AddSecondary(BeamOut->TargetExcitation(),BeamOut->ReactionPosition(),true);

//  	  aParticleChange.DumpInfo();
// 	  getc(stdin);
//  	  aParticleChange.GetSecondary(0)->GetDynamicParticle()->DumpInfo();
//  	  aParticleChange.GetSecondary(1)->GetDynamicParticle()->DumpInfo();
//  	  getc(stdin);
	}
      else
	{
	  aParticleChange.SetNumberOfSecondaries(1);
	  aParticleChange.AddSecondary(BeamOut->ReactionProduct(),BeamOut->ReactionPosition(),true);
	}


    }
   
  return &aParticleChange;
}

G4double Reaction::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
                             G4double,
                             G4ForceCondition* condition
                            )
{


  reaction_here=false;
  *condition=NotForced;

  if(BeamOut->ReactionOn())
    {
      G4String name=aTrack.GetVolume()->GetLogicalVolume()->GetName();
      G4UserLimits* pUserLimits = aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();
      if(name=="target_log")
	{
	  G4double ZReaction=pUserLimits->GetUserMinRange(aTrack);
	  G4double ZCurrent=aTrack.GetPosition().getZ();
	  G4double Z=ZReaction-ZCurrent;
	  if(Z<0)
	    {
// 	      G4cout<<" Past the reaction point"<<G4endl;
// 	      G4cout<<" Volume "<<name<<G4endl;
// 	      G4cout<<" Z[mm]: reaction "<<ZReaction/mm<<" current "<<ZCurrent/mm<<" DZ "<<Z/mm<<G4endl;
	      return DBL_MAX;
	    }
	  if(Z>eps)
	    {
	      G4ThreeVector dir=aTrack.GetDynamicParticle()->GetMomentumDirection();
	      
	      dir*=(ZReaction-ZCurrent);
// 	      G4cout<<" Before the reaction point"<<G4endl;
// 	      G4cout<<" Volume "<<name<<G4endl;
// 	      G4cout<<" Z[mm]: reaction "<<ZReaction/mm<<" current "<<ZCurrent/mm<<" DZ "<<Z/mm<<G4endl;
	      return dir.mag();
	    }
       if(Z<eps)
	 {
	   
	   reaction_here=true;
	   return 0.;
	 }
     }
    }
  return DBL_MAX;
}

