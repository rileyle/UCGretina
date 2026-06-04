#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#include "G4Threading.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4UImanager.hh"

#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "Incoming_Beam.hh"
#include "Incoming_Beam_Messenger.hh"
#include "Outgoing_Beam.hh"
#include "Outgoing_Beam_Messenger.hh"

#ifdef G4VIS_USE
#include "VisManager.hh"
#endif

#include "Git_Hash.hh"

#include "G4Timer.hh"
G4Timer Timer;
G4Timer Timerintern;

int main(int argc,char** argv) 
{
  
  // Construct the default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  G4cout << "Git commit: " << GIT_HASH << G4endl;
  G4cout << "Git branch: " << GIT_BRANCH << G4endl;
  
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

  runManager->SetUserInitialization(new ActionInitialization(detector, BeamIn, BeamOut, /*enableStepping=*/false));

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

  return 0;
}
