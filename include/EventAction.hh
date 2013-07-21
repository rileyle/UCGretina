#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "TrackerIonSD.hh"
#include "TrackerGammaSD.hh"
#include "G4Event.hh"
#include "EventInformation.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

#include "GEB.hh"

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
    void writeDecomp(long long int, G4int, G4int*, G4int*, G4double*, G4double*, G4double*, G4double*);
    void writeSim(long long int, EventInformation*);
    void openCrmatFile();
    void closeCrmatFile();
    void SetCrmatFile(G4String);
    void SetGretinaCoords();
    void SetPackRes(G4double res) { packingRes = res; }
    void SetS800KE(G4double ke) { S800KE = ke; }
    G4double GetS800KE() { return S800KE; }
    const G4Event* GetEvent() { return evt; }
    void SetPrint(){ print = true; }

    void SetInBeam(G4bool flag){fisInBeam = flag;}
    G4bool InBeam(){return fisInBeam;}

    void SetNTotalevents(G4int n){NTotalEvents = n;}
    G4int GetNTotalevents(){return NTotalEvents;}
    void SetEveryNEvents(G4int n){everyNevents = n;}
    G4int GetEveryNEvents(){return everyNevents;}

  private:
    G4int ionCollectionID;
    G4int gammaCollectionID;
    G4String outFileName;
    std::ofstream evfile;
    G4bool evOut;
    G4String crmatFileName;
    G4int crmatFile;
    G4String mode2FileName;
    G4int mode2file;
    G4bool mode2Out;  
    G4float crmat[MAXDETPOS][MAXCRYSTALNO][4][4];
    G4bool gretinaCoords;
    G4double packingRes;
    G4double S800KE;
    const G4Event* evt;
    G4bool print;
    G4bool fisInBeam;
    G4int NTotalEvents;

    G4int timerCount;
    G4int everyNevents;
    G4double eventsPerSecond;
};

#endif //EVENTACTION_H
