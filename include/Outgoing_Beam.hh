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

#include "AngularDistribution.hh"
#include "Randomize.hh"
#include "Incoming_Beam.hh"
#include "Charge_State.hh"



// this is?
#define  eps 0.00001

class Outgoing_Beam
{
public:
  Outgoing_Beam();
  ~Outgoing_Beam();

  void Report();
  void defaultIncomingIon(Incoming_Beam *);
  void setDecayProperties(void);
  void setDA(G4int);
  void setDZ(G4int);
  void setEx(G4double);
  void setLvlSchemeFile(G4String name){lvlSchemeFileName = name;}
  void openLvlSchemeFile();
  void closeLvlSchemeFile();
  void setTarEx(G4double);
  void setTFrac(G4double);
  void settau(G4double);
  void setbeta(G4double);   // all defunct, see .cc file
  void setDoppZ(G4double);  // 
  void setDoppY(G4double);  // 
  void setDoppX(G4double);  // 
  void ScanInitialConditions(const G4Track &);
  void SetReactionOn(){ reacted=true;};
  void SetReactionOff(){reacted=false;}
  void SetSource(){source=true;}
  void SetAlphaTarget(G4double a){alpha=a;Calc_pmax();}         // appear defunct
  void SetThetaMaxTarget(G4double t){theta_max=t;Calc_pmax();}  //
  void SetThetaSigmaA(G4double s){sigma_a=s;} //LR // TB
  void SetThetaSigmaB(G4double s){sigma_b=s;} //LR // TB

  G4ParticleTable* particleTable; 
  G4DynamicParticle* ReactionProduct();
  G4DynamicParticle* ProjectileGS();
  G4DynamicParticle* TargetExcitation();
  G4ThreeVector ReactionPosition();
  G4double GetBetaDopp(){return betaDopp;}
  G4double GetGammaDopp(){return 1./sqrt(1.-betaDopp*betaDopp);}
  G4ThreeVector GetPosDopp(){return posDopp;}
  G4double getTime(){return tau;}
  G4bool   ReactionOn(){return reacted;}
  G4bool   Source(){return source;}
  G4double GetAlpha(){return alpha;}
  G4double GetTheta(){return theta_max;}
  G4double GetThetaSigmaA(){return sigma_a;} //TB
  G4double GetThetaSigmaB(){return sigma_b;} //TB
  G4double getTFrac(){return TFrac;};
  G4int    GetReactionFlag(){return ReactionFlag;}
  void     SetReactionFlag(G4int f){ReactionFlag=f;}
  void     SetNQ(G4int n){NQ=n;SetUpChargeStates();}
  void     SetUpChargeStates();
  void     SelectQ(G4int q){SQ=q;G4cout<<" Charge state "<<SQ<<" selected for setup"<<G4endl;}
  void     SetQCharge(G4int);
  void     SetQUnReactedFraction(G4double);
  void     SetQReactedFraction(G4double);
  void     SetQKEu(G4double);
  void     SetQKE(G4double);
  void     CalcQR();
  void     CalcQUR();
  G4int    GetIndex(){return Index;}
  G4double GetURsetKE();
  G4double GetRsetKE();
 	void SetCoeff(int,double);  // TB 
    void SetTargetCoeff(int, double); // TB 
private:
  G4int Ain;
  G4int Zin;
  G4ThreeVector dirIn;
  G4ThreeVector posIn;
  G4ThreeVector posOut;
  G4ThreeVector pIn;
  G4int  ReactionFlag;
 

  G4double      tauIn;
  G4double      KEIn;

  G4int DZ;
  G4int DA;

  G4double Ex,TarEx,TFrac;
  G4String lvlSchemeFileName;
  std::ifstream lvlSchemeFile;
  G4int    Nlevels;
  G4double relPop[1000];
  G4double levelEnergy[1000];
  G4double tau;
  G4double betaDopp;
  G4int    NQ,SQ;
  vector<Charge_State*> Q;
  G4double  QR[1000],QUR[1000];
  G4int     QRI[1000],QURI[1000],Index;

  G4ThreeVector posDopp;
  G4ParticleDefinition* ion;
  G4ParticleDefinition* iongs;
 
  static G4Decay decay;
  G4bool  reacted;
  G4bool  source;

  Incoming_Beam *beamIn;
  // polynomial coefficients for momentum change for the reacted beam
  G4double alpha;
  G4double sigma_a; //TB
  G4double sigma_b; //TB
  G4double theta_max;
  G4double pmax;
  G4double twopi;

  // TB angular distribution coefficients for a coulex angular distribution
  //	double ai[3];
  AngularDistribution theAngularDistribution;

  G4ThreeVector GetOutgoingMomentum();
  
  // TB angular distribution also for the target.
  //  double targetai[3];
  G4ThreeVector TargetAngularDistribution();
  AngularDistribution theTargetAngularDistribution;

  G4double GetDTheta();
  void Calc_pmax();
};

#endif


           
