
#ifndef TrackerGammaHit_h
#define TrackerGammaHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
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
      void SetDetNumb(G4int num) {detNumb=num;};
      void SetSegNumb(G4int num) {segNumb=num;};
      void SetEdep     (G4double de)      { edep = de; };
      void SetTotalEnergy(G4double te)      { etotal = te; };
      void SetPos      (G4ThreeVector xyz){ pos = xyz; };
      
      G4int GetTrackID()        { return trackID; };
      G4String GetParticleID()  { return particleID; };
      G4String GetProcess()     { return process; };
      G4int  GetDetNumb()       { return detNumb; };
      G4int  GetSegNumb()       { return segNumb; };
      G4double GetEdep()        { return edep; };
      G4double GetTotalEnergy() { return etotal; };
      G4ThreeVector GetPos()    { return pos; };
      
  private:
  
      G4int         trackID;
      G4String      particleID;
      G4String      process;
      G4int         detNumb;
      G4int         segNumb;
      G4double      edep;
      G4double      etotal;
      G4ThreeVector pos;
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
