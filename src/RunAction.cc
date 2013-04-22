#include "RunAction.hh"

#include "G4Timer.hh"
extern G4Timer Timer;

RunAction::RunAction(DetectorConstruction* detector, Outgoing_Beam* BO,EventAction* ev): myDetector(detector), BeamOut(BO),evaction(ev)
{
  
}


RunAction::~RunAction()
{

}

void RunAction::BeginOfRunAction(const G4Run* run)
{

  G4cout<<" Begining of run "<<G4endl;

  evaction->SetNTotalevents(run->GetNumberOfEventToBeProcessed());

  if(BeamOut->ReactionOn()) {
    G4cout<<" Simulating " << run->GetNumberOfEventToBeProcessed()
	  << " scattered ions "<<G4endl;
    evaction->SetInBeam(true);
  } else {
    G4cout<<" Simulating "  << run->GetNumberOfEventToBeProcessed()
	  << " source events "<<G4endl;
    evaction->SetInBeam(false);
  }

  Timer.Start();
}


 
void RunAction::EndOfRunAction(const G4Run*)
{
  if(evaction->EvOut())
    evaction->closeEvfile();
  if(evaction->Mode2Out())
    evaction->closeMode2file();

  Timer.Stop();

  G4cout << "                                                     " << G4endl;

  G4double time, hours, minutes, seconds;

  G4cout << "Real time: ";
  time = Timer.GetRealElapsed();
  hours = floor(time/3600.0);
  if(hours>0){
    G4cout << std::setprecision(0) << std::setw(2) 
	   << hours << ":";
    G4cout << std::setfill('0');
  } else {
    G4cout << std::setfill(' ');
  }
  minutes = floor(fmod(time,3600.0)/60.0);
  if(minutes>0){
    G4cout << std::setprecision(0) << std::setw(2) << minutes << ":";
    G4cout << std::setfill('0');
  } else {
    G4cout << std::setfill(' ');
  }
  seconds = fmod(time,60.0);
  if(seconds>0)
    G4cout << std::setprecision(2) << std::setw(4) << seconds;
  G4cout << std::setfill(' ');

  G4cout << "   System time: ";
  time = Timer.GetSystemElapsed();
  hours = floor(time/3600.0);
  if(hours>0){
    G4cout << std::setprecision(0) << std::setw(2) 
	   << hours << ":";
    G4cout << std::setfill('0');
  } else {
    G4cout << std::setfill(' ');
  }
  minutes = floor(fmod(time,3600.0)/60.0);
  if(minutes>0){
    G4cout << std::setprecision(0) << std::setw(2) << minutes << ":";
    G4cout << std::setfill('0');
  } else {
    G4cout << std::setfill(' ');
  }
  seconds = fmod(time,60.0);
  if(seconds>0)
    G4cout << std::setprecision(2) << std::setw(4) << seconds;
  G4cout << std::setfill(' ');

  G4cout << "   User time: ";
  time = Timer.GetUserElapsed();
  hours = floor(time/3600.0);
  if(hours>0){
    G4cout << std::setprecision(0) << std::setw(2) 
	   << hours << ":";
    G4cout << std::setfill('0');
  } else {
    G4cout << std::setfill(' ');
  }
  minutes = floor(fmod(time,3600.0)/60.0);
  if(minutes>0){
    G4cout << std::setprecision(0) << std::setw(2) << minutes << ":";
    G4cout << std::setfill('0');
  } else {
    G4cout << std::setfill(' ');
  }
  seconds = fmod(time,60.0);
  if(seconds>0)
    G4cout << std::setprecision(2) << std::setw(5) << seconds;
  G4cout << std::setfill(' ');

  G4cout << "   "
	 << evaction->GetNTotalevents()/Timer.GetRealElapsed()
	 << " events/s" << G4endl;
 
}

