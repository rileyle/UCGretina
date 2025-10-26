#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "TrackerIonSD.hh"
#include "TrackerGammaSD.hh"
#include "G4Event.hh"
#include "PrimaryVertexInformation.hh"
#include "AngularDistribution.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

#include "GEB.hh"
#include "cache.hh"

#define MAXDETPOS 31   // for crmat - hole indices offset+1
#define MAXCRYSTALNO 4 // for crmat

#include <fstream>
#include <string>
#include <fcntl.h>

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    void SetOutFile(G4String);
    void SetOutDetsOnly() { outDetsOnly = true; }
    G4String GetOutFileName(){return outFileName;}
    G4bool EvOut(){return evOut;}
    void openEvfile();
    void closeEvfile();
    void SetMode2File(G4String);
    G4String GetMode2FileName(){return mode2FileName;}
    G4bool Mode2Out(){return mode2Out;}
    void openMode2file();
    void closeMode2file();
    void writeGEBHeader(GEBDATA*);
    void writeS800(long long int, G4double, G4double, G4double, G4double);
    void writeDecomp(long long int, G4int, G4int*, G4int*, G4int*, G4double*,
                     G4double*, G4double*, G4double*, G4double*, G4double*);
    void writeSim(long long int, PrimaryVertexInformation*);
    void writeCache(TrackerIonHitsCollection*);
    void openCacheOutputFile(G4String);
    void closeCacheOutputFile();
    G4bool CacheOut(){return cacheOut;}
#ifdef CACHETEXT
    std::ofstream& getCacheOutputFile(){return cacheOutputFile;}
#else
    std::FILE* getCacheOutputFile(){return cacheOutputFile;}
#endif
    void openCacheInputFile(G4String);
    void closeCacheInputFile();
    G4bool CacheIn(){return cacheIn;}
#ifdef CACHETEXT  
    std::ifstream& getCacheInputFile(){return cacheInputFile;}
#else
  std::FILE* getCacheInputFile(){return cacheInputFile;}
#endif
    void SetCacheHalfLife(G4double hl) {cacheHalfLife=hl;}
    void SetCacheGammaEnergy(G4double e) {cacheGammaEnergy=e;}
    void SetCacheZOffset(G4double z) {cacheZOffset=z;}
    G4double GetCacheHalfLife() { return cacheHalfLife; }
    G4double GetCacheGammaEnergy() { return cacheGammaEnergy; }
    G4double GetCacheZOffset() { return cacheZOffset; }
    void SetCacheAngDistA0(G4double);
    void SetCacheAngDistA2(G4double);
    void SetCacheAngDistA4(G4double);
    AngularDistribution* GetCacheAngDist(){return cacheAngDist;}
    void openCrmatFile();
    void closeCrmatFile();
    void SetCrmatFile(G4String);
    void SetCrystalXforms();
    void SetGretinaCoords();
    void SetPackRes(G4double res) { packingRes = res; }
    void SetS800KE(G4double ke) { S800KE = ke; }
    G4double GetS800KE() { return S800KE; }
    void SetAllS800() { allS800 = true; }
    G4bool AllS800() { return allS800; }
    const G4Event* GetEvent() { return evt; }
    void SetTimeSort(){ timeSort = true; }
    void SetPrint(){ print = true; }

    void SetInBeam(G4bool flag){fisInBeam = flag;}
    G4bool InBeam(){return fisInBeam;}

    void SetNTotalevents(G4int n){NTotalEvents = n;}
    G4int GetNTotalevents(){return NTotalEvents;}
    void SetEveryNEvents(G4int n){everyNevents = n;}
    G4int GetEveryNEvents(){return everyNevents;}

    void SetPosRes(G4double res){posRes = res;}
    void SetThreshE(G4double e){threshE = e;}
    void SetThreshDE(G4double de){threshDE = de;}
  
  private:
    G4int ionCollectionID;
    G4int gammaCollectionID;
    G4String outFileName;
    G4bool   outDetsOnly;
    std::ofstream evfile;
    G4bool evOut;
    G4String crmatFileName;
    G4int crmatFile;
    G4String cacheOutputFileName;
#ifdef CACHETEXT
    std::ofstream cacheOutputFile;
#else
  std::FILE* cacheOutputFile;
#endif
    G4bool cacheOut;
    G4String cacheInputFileName;
#ifdef CACHETEXT
    std::ifstream cacheInputFile;
#else
  std::FILE* cacheInputFile;
#endif
    G4bool cacheIn;
    G4double cacheGammaEnergy;
    G4double cacheZOffset;
    G4double cacheHalfLife;
    AngularDistribution* cacheAngDist;
    G4String mode2FileName;
    G4int mode2file;
    G4bool mode2Out;  
    G4float crmat[MAXDETPOS][MAXCRYSTALNO][4][4];
    G4bool gretinaCoords;
    G4bool crystalXforms;
    G4double hitTolerance;
    G4double packingRes;
    G4double S800KE;
    G4bool allS800;
    const G4Event* evt;
    G4bool print;
    G4bool fisInBeam;
    G4bool timeSort;
    G4int NTotalEvents;
    G4double posRes;
    G4double threshE;
    G4double threshDE;
  
    G4int timerCount;
    G4int everyNevents;
    G4double eventsPerSecond;
};

#endif //EVENTACTION_H
