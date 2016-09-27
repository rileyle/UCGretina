#ifndef Outgoing_Beam_h
#define Outgoing_Beam_h 1

#include <cmath>
#include <vector>
using namespace std;

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RandomDirection.hh"
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecaymessenger.hh"
#include "NuclearLevelData.hh"
#include "LevelManager.hh"

#include "Randomize.hh"
#include "Incoming_Beam.hh"

// this is?
#define  eps 0.00001

class Outgoing_Beam
{
public:
  Outgoing_Beam();
  ~Outgoing_Beam();

  void Report();
  void Update(){ beamIn->Report(); setDecayProperties(); Report(); }
  void defaultIncomingIon(Incoming_Beam *);
  void setDecayProperties(void);
  void setDA(G4int);
  void setDZ(G4int);
  void setEx(G4double);
  void setProjectileExcitation(){ targetExcitation = false; }
  void setTargetExcitation(){ targetExcitation = true; }
  G4bool TargetExcitation(){ return targetExcitation; }
  void setLvlDataFile(G4String name){ lvlDataFileName = name; }
  void setTarA(G4int);
  void setTarZ(G4int);
  void setTarEx(G4double);
  void setTFrac(G4double);
  void ScanInitialConditions(const G4Track &);
  void SetReactionOn(){ reacted=true;};
  void SetReactionOff(){reacted=false;}
  void SetSource(){source=true;}
  void SetThetaMin(G4double t){theta_min=t;}
  void SetThetaMax(G4double t){theta_max=t;}
  void SetThetaSigmaA(G4double sig){sigma_a=sig;} //LR // TB
  void SetThetaSigmaB(G4double sig){sigma_b=sig;} //LR // TB

  G4ParticleTable* particleTable; 
  G4DynamicParticle* ReactionProduct();
  G4ThreeVector ReactionPosition();
  G4int getTarA(){return TarA;}
  G4int getTarZ(){return TarZ;}
  G4bool   ReactionOn(){return reacted;}
  G4bool   Source(){return source;}
  G4double GetThetaMax(){return theta_max;}
  G4double GetThetaMin(){return theta_min;}
  G4double GetThetaSigmaA(){return sigma_a;} //TB
  G4double GetThetaSigmaB(){return sigma_b;} //TB
  void setXsectFile(G4String);
  G4double getTFrac(){return TFrac;};
  G4int    GetReactionFlag(){return ReactionFlag;}
  void     SetReactionFlag(G4int f){ReactionFlag=f;}
  G4int    AboveThreshold(){return ThresholdFlag;}

private:
  G4int Ain;
  G4int Zin;
  G4ThreeVector dirIn;
  G4ThreeVector posIn;
  G4ThreeVector posOut;
  G4ThreeVector pIn;
  G4int  ReactionFlag;
  G4int  ThresholdFlag;
  G4double      KEIn;

  G4int DZ;
  G4int DA;
  G4int TarA;
  G4int TarZ;
  G4double m1;
  G4double m2;
  G4double m3;
  G4double m4;
  G4double ET;
  G4double p1;
  G4double sin2theta3_max;

  G4double Ex,TarEx,TFrac;
  G4String lvlDataFileName;
  std::ifstream lvlDataFile;

  G4ParticleDefinition* beam;
  G4ParticleDefinition* ion;
  G4ParticleDefinition* ionGS;
  G4ParticleDefinition* tarIn;
  G4ParticleDefinition* tarOut;
  G4ParticleDefinition* tarOutGS;
  
  static G4RadioactiveDecay* decay;
  static G4RadioactiveDecaymessenger* decayMessenger;
  G4bool  reacted;
  G4bool  source;

  G4bool targetExcitation;

  Incoming_Beam *beamIn;

  G4double sigma_a;
  G4double sigma_b;
  G4String xsectFileName;
  G4double theta_min;
  G4double theta_max;
  G4double theta_bin;
  G4int    Nxsect;
  G4double Xsect[1000];
  G4double twopi;

  G4ThreeVector GetOutgoingMomentum();
  
  G4double GetDTheta();

};

#endif


           
