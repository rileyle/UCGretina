#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
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
  void SetSource(){source=true;inbeam=false;background=false;}
  void SetInBeam(){source=false;inbeam=true;background=false;}
  void SetBackground(){source=false;inbeam=false;background=true;}
  void SetSourceX(G4double x){sourcePosition.setX(x);}
  void SetSourceY(G4double y){sourcePosition.setY(y);}
  void SetSourceZ(G4double z){sourcePosition.setZ(z);}
  void SetSourceR(G4double r){sourceRadius=r;}
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
  void SetSourceRa226();
  void SetSourceRa226Peaks();
  void SetSourceAm241();
  void SetSourcePhotopeaks();
  void SetSourceAu();
  void SetSourceSimple();
  void SetSourceWhite();
  void SetSourceBG();
  void SetSourceBGWhite();
  void SetSourceMuon();
  //  void ReactionOn(){BeamOut->SetReactionOn();fracOn=false;}   //REMOVE
  //  void ReactionOff(){BeamOut->SetReactionOff();fracOn=false;} //REMOVE
  //  void SetFraction(G4double f){fracOn=true;frac=f;} //REMOVE
  void SetSourceEnergy(G4double);
  void SetWhiteSourceLowE(G4double energy){sourceWhiteLoE = energy;}
  void SetWhiteSourceHighE(G4double energy){sourceWhiteHiE = energy;}
  void SetSourceCollAngle(G4double ang){collimationAngle = ang; isCollimated = true;}
  void SetSourceCollDirection(G4ThreeVector dir){collimationDirection = dir; isCollimated = true;}
  G4double GetSourceEnergy();
  G4bool GetInBeam(){return inbeam;}

private:


  G4int n_particle;
  G4ParticleGun* particleGun;
  G4ParticleTable* particleTable;
  G4IonTable* ionTable;
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
  G4bool source, inbeam, background, isCollimated; 
  G4String sourceType;
  G4ThreeVector sourcePosition;
  G4double sourceRadius;
  vector<SourceData*> TheSource;
  G4double sourceBranchingSum;
  G4double sourceWhiteLoE;
  G4double sourceWhiteHiE;
  G4double bgSphereRmin;
  G4double bgSphereRmax;
  G4double collimationAngle;
  G4ThreeVector collimationDirection;
};


#endif
