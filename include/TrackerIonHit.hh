
#ifndef TrackerIonHit_h
#define TrackerIonHit_h 1

#define GUN_FLAG           1
#define TARGET_FACE_FLAG   2
#define TARGET_BACK_FLAG   3
#define REACTION_FLAG      4
#define DECAY_FLAG         5
#define MAX_FLAGS          6

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iomanip>



class TrackerIonHit : public G4VHit
{
  public:

      TrackerIonHit();
     ~TrackerIonHit();
      TrackerIonHit(const TrackerIonHit&);
      const TrackerIonHit& operator=(const TrackerIonHit&);
      G4int operator==(const TrackerIonHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
      void SetTrackID  (G4int track)      { trackID = track; };
      void SetParticleID (G4String particle)      { particleID = particle; }; 
      void SetBeta(G4double b){beta=b;};
      void SetTheta(G4double th){theta=th;};
      void SetPhi(G4double ph){phi=ph;};
      void SetTime(G4double pt){time=pt;};
      void SetLabTime(G4double lt){labtime=lt;};
      void SetGlobTime(G4double gt){globaltime=gt;};
      void SetKE(G4double k){KE=k;};
      void SetEdep(G4double E){Edep=E;};
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };
      void SetGunFlag(){flag=GUN_FLAG;}
      void SetDecayFlag(){flag=DECAY_FLAG;}
      void SetReactionFlag(){flag=REACTION_FLAG;}
      void SetTargetFaceFlag(){flag=TARGET_FACE_FLAG;}
      void SetTargetBackFlag(){flag=TARGET_BACK_FLAG;}
      void SetVolName (G4String vn){volname=vn;};
      void SetLength(G4double len){length=len;}
   
      G4int GetTrackID()    { return trackID; };
      G4String GetParticleID() {return particleID;};
      G4double GetBeta(){return beta;};
      G4double GetTheta(){return theta;};
      G4double GetPhi(){return phi;};
      G4double GetTime(){return time;};
      G4double GetLabTime(){return labtime;};
      G4double GetGlobTime(){return globaltime;};
      G4double GetKE(){return KE;};
      G4double GetEdep(){return Edep;};
      G4ThreeVector GetPos(){ return pos; };
      G4String GetVolName() {return volname;};
      G4int    GetFlag(){return flag;}
      G4double GetLength(){return length;}
  private:
  
      G4int         trackID;
      G4String      particleID;
      G4double      beta;
      G4double      theta;
      G4double      phi;
      G4double      time;
      G4double      labtime;
      G4double      globaltime;
      G4double      KE;
      G4double      Edep;
      G4ThreeVector pos;
      G4String      volname;
      G4int         flag;
      G4double      length;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TrackerIonHit> TrackerIonHitsCollection;

extern G4Allocator<TrackerIonHit> TrackerIonHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TrackerIonHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) TrackerIonHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TrackerIonHit::operator delete(void *aHit)
{
  TrackerIonHitAllocator.FreeSingle((TrackerIonHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
