#include "RunAction.hh"
#include "cache.hh"
#include "G4Timer.hh"
extern G4Timer Timer;

RunAction::RunAction(DetectorConstruction* detector, Incoming_Beam* BI,EventAction* ev): myDetector(detector), BeamIn(BI), evaction(ev)
{

}


RunAction::~RunAction()
{

}

void RunAction::BeginOfRunAction(const G4Run* run)
{

  G4cout<<" Beginning of run "<<G4endl;

  evaction->SetNTotalevents(run->GetNumberOfEventToBeProcessed());
  if(run->GetNumberOfEventToBeProcessed() > 1000000)
    evaction->SetEveryNEvents(10000);
  else if(run->GetNumberOfEventToBeProcessed() > 1000)
    evaction->SetEveryNEvents(1000);
  else if(run->GetNumberOfEventToBeProcessed() > 100)
    evaction->SetEveryNEvents(100);
  else
    evaction->SetEveryNEvents(1);

  G4cout << " Simulating " << run->GetNumberOfEventToBeProcessed()
	 << " events." << G4endl;

  if(BeamIn->getKE()>0)
    evaction->SetInBeam(true);
  
  if(evaction->EvOut())
    G4cout << " Writing ASCII output to " 
	   << evaction->GetOutFileName() << G4endl;
  if(evaction->Mode2Out())
    G4cout << " Writing Mode 2 output to " 
	   << evaction->GetMode2FileName() << G4endl;
  if(evaction->CacheOut()){
#ifdef CACHETEXT
    std::ofstream& cacheOutputFile = evaction->getCacheOutputFile(); 
    G4int Nevents = run->GetNumberOfEventToBeProcessed();
    cacheOutputFile << Nevents << G4endl;
#else
    std::FILE* cacheOutputFile = evaction->getCacheOutputFile(); 
    G4int Nevents = run->GetNumberOfEventToBeProcessed();
    fwrite(&Nevents, sizeof(G4int), 1, cacheOutputFile);
#endif
  }
  if(evaction->CacheIn()){
#ifdef CACHETEXT
    std::ifstream& cacheInputFile = evaction->getCacheInputFile();
    if(cacheInputFile.eof()){
      G4cerr << "Error reading cache file header." << G4endl;
      exit(EXIT_FAILURE);
    }
    G4int Nevents;
    cacheInputFile >> Nevents;
#else
    std::FILE* cacheInputFile = evaction->getCacheInputFile();
    G4int Nevents;
    int size = fread(&Nevents, sizeof(G4int), 1, cacheInputFile);
    if(size != 1) {
      G4cerr << "Error reading cache file header." << G4endl;
      exit(EXIT_FAILURE);
    }
#endif
    G4cout << Nevents << " events in cache file." << G4endl;
    if (Nevents < run->GetNumberOfEventToBeProcessed()){
      G4cerr << "Error: There are only " << Nevents
	     << " events in the cache file and the user has requested "
	     << run->GetNumberOfEventToBeProcessed() << G4endl;
      exit(EXIT_FAILURE);
    }
  }
  Timer.Start();
}


 
void RunAction::EndOfRunAction(const G4Run*)
{
  if(evaction->EvOut())
    evaction->closeEvfile();
  if(evaction->Mode2Out())
    evaction->closeMode2file();
  if(evaction->CacheOut())
    evaction->closeCacheOutputFile();
  if(evaction->CacheIn())
    evaction->closeCacheInputFile();

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
    G4cout << std::setprecision(2) << std::setw(5) << seconds;
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
    G4cout << std::setprecision(2) << std::setw(5) << seconds;
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

