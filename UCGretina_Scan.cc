#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorAction_Messenger.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "EventAction_Messenger.hh"
#include "RunAction.hh"
#include "Incoming_Beam.hh"
#include "Incoming_Beam_Messenger.hh"
#include "Outgoing_Beam.hh"
#include "Outgoing_Beam_Messenger.hh"

#ifdef G4VIS_USE
#include "VisManager.hh"
#endif

#include "G4Timer.hh"
G4Timer Timer;
G4Timer Timerintern;

int main(int argc,char** argv) 
{
  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  cout << "Instantiating DetectorConstruction ..." << endl;
  // set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

  PhysicsList *physicsList = new PhysicsList(detector);
  runManager->SetUserInitialization(physicsList);

  cout << "... Done" << endl;

  // Construct incoming and outgoing beams
  Incoming_Beam* BeamIn = new Incoming_Beam();
  Incoming_Beam_Messenger* IncomingBeamMessenger = new Incoming_Beam_Messenger(BeamIn);

  Outgoing_Beam* BeamOut = new Outgoing_Beam();
  BeamOut->defaultIncomingIon(BeamIn);
  physicsList->SetOutgoingBeam(BeamOut);
  Outgoing_Beam_Messenger* OutgoingBeamMessenger = new Outgoing_Beam_Messenger(BeamOut);

  // set mandatory user action class
  EventAction* eventAction = new EventAction();
  EventAction_Messenger* eventActionMessenger = new EventAction_Messenger(eventAction);
  runManager->SetUserAction(eventAction);

  PrimaryGeneratorAction* generatorAction = new PrimaryGeneratorAction(detector,BeamIn,BeamOut);
  PrimaryGeneratorAction_Messenger* generatorActionMessenger = new PrimaryGeneratorAction_Messenger(generatorAction);

  runManager->SetUserAction(generatorAction);
  RunAction* runAction = new RunAction(detector,BeamIn,eventAction);
  runManager->SetUserAction(runAction);

  TrackingAction* trackingAction = new TrackingAction(eventAction);
  runManager->SetUserAction(trackingAction);

  //  User Stepping Action only needed for in-beam simulations.
  //  SteppingAction* steppingAction = new SteppingAction();
  //  runManager->SetUserAction(steppingAction);

  G4UIsession* session=0;

#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager=0;
#endif

  if (argc==1)   // Define UI session for interactive mode.
    {

#ifdef G4VIS_USE
      // visualization manager
      cout << "Starting visualization...";
      visManager = new VisManager; 
      visManager->Initialize();
      cout << "Done!" << endl;
#endif

// G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);      
#else
      session = new G4UIterminal();
#endif

    }

  // Initialize G4 kernel
  // cout << "*** Initializing runManager" << endl;
  // //  runManager->SetVerboseLevel(2);
  // runManager->Initialize();
  // cout << "*** Initialized runManager" << endl;

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (session)   // Define UI session for interactive mode.
    {
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    { 
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  // job termination
  if(argc==1){
#ifdef G4VIS_USE
  delete visManager;
#endif
  }

  delete runManager;

  delete BeamIn;

  delete IncomingBeamMessenger;

  delete BeamOut;

  delete OutgoingBeamMessenger;

  delete eventActionMessenger;

  delete generatorActionMessenger;

  return 0;
}
