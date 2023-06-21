
#ifndef Reaction_h
#define Reaction_h 1

#include <unordered_map>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VProcess.hh"
#include "Outgoing_Beam.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
#include "G4UserLimits.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4NuclearPolarizationStore.hh"
#include "G4Fragment.hh"
#include "G4Clebsch.hh"

#define  eps 0.00001

class Reaction_Messenger;

class Reaction : public G4VProcess 
{
  public:     
   G4bool reaction_here;
   G4bool decayed_at_rest;
   G4bool target_reaction;
   G4bool ground_state;

     Reaction(Outgoing_Beam*, const G4String& processName ="Reaction" );

     virtual ~Reaction();

 
 
    virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    );

     virtual G4VParticleChange* PostStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    );
			    
     //  no operation in  AtRestGPIL
     virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& ,
			     G4ForceCondition* 
			    ){ return DBL_MAX; };
			    
     //  no operation in  AtRestDoIt      
     virtual G4VParticleChange* AtRestDoIt(
			     const G4Track& ,
			     const G4Step&
			    );
			     //			    ){ return NULL; };

     //  no operation in  AlongStepGPIL
     virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track&,
			     G4double  ,
			     G4double  ,
			     G4double& ,
                             G4GPILSelection*
			    ){ return -1.0; };

     //  no operation in  AlongStepDoIt
     virtual G4VParticleChange* AlongStepDoIt(
			     const G4Track& ,
			     const G4Step& 
			    ) {return NULL;};

  void SetPopulation(G4int,G4double);

  private:
  
  // hide assignment operator as private 
  Reaction& operator=(const Reaction&){return *this;};

  Outgoing_Beam* BeamOut;
  Reaction_Messenger *theMessenger;
  std::unordered_map<int,double> substates;
};

#endif

