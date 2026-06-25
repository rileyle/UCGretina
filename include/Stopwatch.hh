#pragma once
#include "G4Timer.hh"
#include <atomic>

struct Stopwatch{
  G4int timerCount = 0;
  G4Timer* Timer = nullptr;
  G4double rate = 0;
  std::atomic<G4int> completedEvents{0};
  
};
