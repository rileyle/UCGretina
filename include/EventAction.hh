#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "TrackerIonSD.hh"
#include "TrackerGammaSD.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include <fstream>
#include <string>
//#include "Results.hh"
class EventAction : public G4UserEventAction
{
  public:
    EventAction();      //LR EventAction(Results*);
   ~EventAction();

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    void SetOutFile(G4String);
    void openEvfile();
    void closeEvfile();

    void openCrmatFile();
    void closeCrmatFile();
    void SetCrmatFile(G4String);
    void SetGretinaCoords();
    void SetPosRes(G4double res) { positionRes = res; };

  private:
    G4int ionCollectionID;
    G4int gammaCollectionID;
    G4String outFileName;
    std::ofstream evfile;
    G4String crmatFileName;
    std::ifstream crmatFile;
    G4double crmat[30][4][4][4];
    G4bool gretinaCoords;
    G4double positionRes;
};

#endif //EVENTACTION_H
