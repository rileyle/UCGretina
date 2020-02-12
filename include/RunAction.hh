#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "DetectorConstruction.hh"
#include "G4Run.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"
#include "Incoming_Beam.hh"
#include "EventAction.hh"

class DetectorConstruction;

class RunAction : public G4UserRunAction
{
  public:
    RunAction(DetectorConstruction*, Incoming_Beam*, EventAction*);
   ~RunAction();

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);

  private:
  DetectorConstruction* myDetector;
  Incoming_Beam* BeamIn;
  EventAction* evaction;

};




#endif

    
