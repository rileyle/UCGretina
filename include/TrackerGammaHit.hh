
#ifndef TrackerGammaHit_h
#define TrackerGammaHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class TrackerGammaHit : public G4VHit
{
  public:

      TrackerGammaHit();
     ~TrackerGammaHit();
      TrackerGammaHit(const TrackerGammaHit&);
      const TrackerGammaHit& operator=(const TrackerGammaHit&);
      G4int operator==(const TrackerGammaHit&) const;

      inline void* operator new(size_t);
      inline void  operator delete(void*);

      void Draw();
      void Print();

  public:
  
      void SetTrackID  (G4int track)         { trackID = track; };
      void SetParticleID (G4String particle) { particleID = particle; };
      void SetProcess (G4String proc)        { process = proc; };
      void SetParentTrackID (G4int parent)   { parentTrackID = parent; };
      void SetCreatorProcess (G4String proc) { creatorProcess = proc; };
      void SetDetNumb(G4int num) {detNumb=num;};
      void SetSegNumb(G4int num) {segNumb=num;};
      void SetEdep     (G4double de)      { edep = de; };
      void SetKE       (G4double e)       { ke = e; };
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };
      void SetPosCrys  (G4ThreeVector xyz){ posCrys = xyz; };
      void SetTrackOrigin(G4ThreeVector xyz){ trackOrigin = xyz; };
      void SetGlobalTime(G4double t)       { globalTime = t; };

      G4int GetTrackID()             { return trackID; };
      G4String GetParticleID()       { return particleID; };
      G4String GetProcess()          { return process; };
      G4int GetParentTrackID()       { return parentTrackID; };
      G4String GetCreatorProcess()   { return creatorProcess; };
      G4int  GetDetNumb()            { return detNumb; };
      G4int  GetSegNumb()            { return segNumb; };
      G4double GetEdep()             { return edep; };
      G4double GetKE()               { return ke; };
      G4ThreeVector GetPos()         { return pos; };
      G4ThreeVector GetPosCrys()     { return posCrys; };
      G4ThreeVector GetTrackOrigin() { return trackOrigin; };
      G4double GetGlobalTime()       { return globalTime; };
      
  private:
  
      G4int         trackID;
      G4String      particleID;
      G4String      process;
      G4int         parentTrackID;
      G4String      creatorProcess;
      G4int         detNumb;
      G4int         segNumb;
      G4double      edep;
      G4double      ke;
      G4ThreeVector pos;
      G4ThreeVector posCrys;
      G4ThreeVector trackOrigin;
      G4double      globalTime;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<TrackerGammaHit> TrackerGammaHitsCollection;

extern G4Allocator<TrackerGammaHit> TrackerGammaHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* TrackerGammaHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) TrackerGammaHitAllocator.MallocSingle();
  return aHit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void TrackerGammaHit::operator delete(void *aHit)
{
  TrackerGammaHitAllocator.FreeSingle((TrackerGammaHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
