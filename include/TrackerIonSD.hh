
#ifndef TrackerIonSD_h
#define TrackerIonSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TrackerIonHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4VTouchable.hh"
#include "G4VProcess.hh"
#include "G4StepStatus.hh"
#include "G4TrackStatus.hh"
#include "G4RunManager.hh"

class TrackerIonSD : public G4VSensitiveDetector
{
  public:
      TrackerIonSD(G4String);
     ~TrackerIonSD();

      void SetPrint(){
	G4cout<<"----> Ion track data set to print at the end of event"<<G4endl;
	print=true;}
      void UnSetPrint(){
	G4cout<<"----> Ion track data set not to print at the end of event"<<G4endl;
	print=false;}

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

 
  private:
      TrackerIonHitsCollection* ionCollection;
      G4bool print;
    
};


#endif
