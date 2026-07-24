#include "RunAction.hh"
#include "cache.hh"
#include "G4Timer.hh"
extern G4Timer Timer;

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif

RunAction::RunAction(DetectorConstruction* detector, Incoming_Beam* BI,EventAction* ev): myDetector(detector), BeamIn(BI), evaction(ev)
{

}


RunAction::~RunAction()
{

}

void RunAction::BeginOfRunAction(const G4Run* run)
{

  if(G4Threading::IsMasterThread()){
    
    G4cout<<" Beginning of run "<<G4endl;

    G4cout << " Simulating " << run->GetNumberOfEventToBeProcessed()
	   << " events." << G4endl;
  
    if(evaction->EvOut())
      G4cout << " Writing ASCII output to " 
	     << evaction->GetOutFileName() << G4endl;
    if(evaction->Mode2Out())
      G4cout << " Writing Mode 2 output to " 
	     << evaction->GetMode2FileName() << G4endl;
    Timer.Start();

  }
  // Write a "rough draft" header to the cache file, to be replaced
  // when we know how many trajectories were written.
  if(evaction->CacheOut()){
    G4int nThreads = 1;
    
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = G4MTRunManager::GetMasterRunManager();
    nThreads = runManager->GetNumberOfThreads();
#endif

    G4int eventsToBeProcessed = float(run->GetNumberOfEventToBeProcessed())/float(nThreads);

#ifdef CACHETEXT
    std::ofstream& cacheOutputFile = evaction->getCacheOutputFile(); 
    //    G4int Nevents = run->GetNumberOfEventToBeProcessed(); 
    cacheOutputFile << eventsToBeProcessed << G4endl;
#else
    std::FILE* cacheOutputFile = evaction->getCacheOutputFile(); 
    //    G4int Nevents = run->GetNumberOfEventToBeProcessed();
    //    fwrite(&Nevents, sizeof(G4int), 1, cacheOutputFile);
    fwrite(&eventsToBeProcessed, sizeof(G4int), 1, cacheOutputFile);
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

    G4int nThreads = 1;
    
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = G4MTRunManager::GetMasterRunManager();
    nThreads = runManager->GetNumberOfThreads();
#endif

    G4int eventsToBeProcessed = float(run->GetNumberOfEventToBeProcessed())/float(nThreads);
    
    if (Nevents < eventsToBeProcessed){
      G4cerr << "Error: There are only " << Nevents
	     << " events in the cache file";
      G4cerr << ", and the user has requested "
	     << eventsToBeProcessed;
#ifdef G4MULTITHREADED
      G4cerr << " per thread";
#endif
      G4cerr << "." << G4endl;
      exit(EXIT_FAILURE);
    }
  }

  evaction->SetNTotalevents(run->GetNumberOfEventToBeProcessed());
  if(run->GetNumberOfEventToBeProcessed() > 1000000)
    evaction->SetEveryNEvents(10000);
  else if(run->GetNumberOfEventToBeProcessed() > 1000)
    evaction->SetEveryNEvents(1000);
  else if(run->GetNumberOfEventToBeProcessed() > 100)
    evaction->SetEveryNEvents(100);
  else
    evaction->SetEveryNEvents(1);

  if(BeamIn->getKE()>0)
    evaction->SetInBeam(true);

  
}


 
void RunAction::EndOfRunAction(const G4Run* run)
{

  if(evaction->CacheOut()){
    evaction->closeCacheOutputFile(); //close it so we can reopen and rewrite header
#ifdef CACHETEXT
    std::fstream cacheOutputFile(evaction->GetCacheOutputFilename(), std::ios::in | std::ios::out);
    cacheOutputFile.seekp(0, std::ios::beg); 
    G4int Nevents = evaction->GetCompletedEvents();
    cacheOutputFile << Nevents << G4endl;
    cacheOutputFile.close();
#else
    std::FILE* cacheOutputFile = std::fopen(evaction->GetCacheOutputFilename(), "r+");
    std::fseek(cacheOutputFile, 0, SEEK_SET);
    G4int Nevents = evaction->GetCompletedEvents();
    fwrite(&Nevents, sizeof(G4int), 1, cacheOutputFile);
    std::fclose(cacheOutputFile);
#endif
  }

  if(evaction->EvOut())
    evaction->closeEvfile();
  if(evaction->Mode2Out())
    evaction->closeMode2file();
  if(evaction->CacheIn())
    evaction->closeCacheInputFile();

  if(G4Threading::IsMasterThread()){
    
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
	 << run->GetNumberOfEventToBeProcessed()/Timer.GetRealElapsed()
	 << " events/s" << G4endl;

  }
  
}

