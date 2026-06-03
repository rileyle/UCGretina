#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "Incoming_Beam.hh"
#include "Outgoing_Beam.hh"

#include "EventAction.hh"
#include "EventAction_Messenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction_Messenger.hh"
#include "RunAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization(DetectorConstruction* detector,
                                           Incoming_Beam* beamIn,
                                           Outgoing_Beam* beamOut,
                                           bool enableStepping)
  : fDetector(detector), fBeamIn(beamIn), fBeamOut(beamOut), fEnableStepping(enableStepping) {}

void ActionInitialization::BuildForMaster() const {
  // Master does not process events; keep it minimal.
  G4cout << "Building Master" << G4endl;
  auto* eventAction = new EventAction();
  SetUserAction(new RunAction(fDetector, fBeamIn, eventAction));
}

void ActionInitialization::Build() const {
  auto* eventAction = new EventAction();
  SetUserAction(eventAction);
  (void)new EventAction_Messenger(eventAction);

  auto* generatorAction = new PrimaryGeneratorAction(fDetector, fBeamIn, fBeamOut);
  SetUserAction(generatorAction);
  (void)new PrimaryGeneratorAction_Messenger(generatorAction);

  SetUserAction(new RunAction(fDetector, fBeamIn, eventAction));
  SetUserAction(new TrackingAction(eventAction));
  if (fEnableStepping) {
    SetUserAction(new SteppingAction());
  }
}
