#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "DetectorConstruction.hh"
#include "Incoming_Beam.hh"
#include "Outgoing_Beam.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include <vector>
#include "SourceData.hh"

using namespace std;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(DetectorConstruction*,Incoming_Beam*,Outgoing_Beam*);
  ~PrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  void SetSource(){source=true;inbeam=false;};
  void SetInBeam(){source=false;inbeam=true;};
  void SetSourceX(G4double x){sourcePosition.setX(x);}
  void SetSourceY(G4double y){sourcePosition.setY(y);}
  void SetSourceZ(G4double z){sourcePosition.setZ(z);}
  void SetSourceOnTargetFace();
  void SetSourceOnTargetBack();
  void SourceReport();
  void SetSourceType(G4String name);
  void SetSourceEu152();
  void SetSourceEu152Peaks();
  void SetSourceCs137();
  void SetSourceCo56();
  void SetSourceCo56Peaks();
  void SetSourceCo60();
  void SetSourcePhotopeaks();
  void SetSourceAu();
  void SetSourceSimple();
  void ReactionOn(){BeamOut->SetReactionOn();fracOn=false;}
  void ReactionOff(){BeamOut->SetReactionOff();fracOn=false;}
  void SetFraction(G4double f){fracOn=true;frac=f;}
  G4double GetSourceEnergy();
  void SetSourceEnergy(G4double);

private:


  G4int n_particle;
  G4ParticleGun* particleGun;
  G4ParticleTable* particleTable;
  G4ParticleDefinition* ion;
  G4ThreeVector  direction;
  G4ThreeVector  position;
  G4double       KE;
  Incoming_Beam* BeamIn;
  Outgoing_Beam* BeamOut;
  DetectorConstruction* myDetector;
  G4double frac;
  G4bool   fracOn;
  // source stuff
  G4bool source,inbeam; 
  G4String sourceType;  //LR
  G4ThreeVector sourcePosition;
  vector<SourceData*> TheSource;
  G4double sourceBranchingSum;
};


#endif


           
