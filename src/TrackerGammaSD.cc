
#include "TrackerGammaSD.hh"
#include "G4RunManager.hh"
#include "DetectorConstruction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4VTouchable.hh"
#include <string.h>
#include <stdio.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerGammaSD::TrackerGammaSD(G4String name)
:G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="gammaCollection");
  print = false;
  trackingMethod = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackerGammaSD::~TrackerGammaSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerGammaSD::Initialize(G4HCofThisEvent*)
{

  gammaCollection = new TrackerGammaHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  annihilationHappened = false; 
  residualEnergy = 0.;  

  G4RunManager* runManager = G4RunManager::GetRunManager();
  if(runManager->GetCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetG4code()->GetParticleName()
     == "gamma")
    primaryGammaID = 0; // Source and background simulations
  else
    primaryGammaID = 2; // In-beam simulations

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TrackerGammaSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

  G4String name = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if(name.substr(0,5)=="geCap")
    {

      G4RunManager* runManager = G4RunManager::GetRunManager();
      DetectorConstruction* theDetector = (DetectorConstruction*)runManager->GetUserDetectorConstruction();

      // Code for trackingMethod case 1 adopted from the 2009 release of the AGATA simulation code.
      G4double edep = 0.;
      G4double additionalEnergy = 0.;
      G4double deltaEnergy = 0.;
      G4String ChiEStato;
      G4String ilProcesso;
      G4String ilGeneratore = " ";
      G4int padre;

      switch(trackingMethod){
      case 1:
	///////////////////////////////////////////////////////////////////
	//// Idea: primary gamma is followed
	//// store the position where an e+ annihilated
	//// accept only secondary gammas coming from that position!!!
	///////////////////////////////////////////////////////////////////

	/// Consider only gammas and positrons
	ChiEStato = aStep->GetTrack()->GetDefinition()->GetParticleName();
	ilProcesso = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	//	G4cout << ChiEStato << " interacted through a " << ilProcesso << G4endl;
	/// positrons: just store the position!!!
	if( ChiEStato == "e+" && ilProcesso == "annihil" ) {
	  positionOfAnnihilation = aStep->GetPostStepPoint()->GetPosition();
	  /// some positron decays in flight, should correct for that ...
	  //	  residualEnergy = aStep->GetTrack()->GetKineticEnergy(); // LR: always zero?
	  // LR: This gives us our escape peaks:
	  residualEnergy = -1.*(aStep->GetPostStepPoint()->GetKineticEnergy() 
				- aStep->GetPreStepPoint()->GetKineticEnergy());
	  annihilationHappened = true;
	  //	  G4cout << " e+ annihilated at " << positionOfAnnihilation/mm << G4endl; 
	  //	  G4cout << " e+ residual energy is " << residualEnergy/keV << G4endl;
	  if( residualEnergy ) {
	    //	    G4cout << ChiEStato << " interacted through a " << ilProcesso << G4endl;
	    //	    G4cout << " e+ residual energy is " << residualEnergy/keV << G4endl;
	    G4int ii = gammaCollection->entries();
	    //	    G4cout << " number of entries is " << ii << G4endl;
	    if( ii ) {
	      //	      G4cout << " Deposited energy was " << (*gammaCollection)[ii-1]->GetEdep()/keV << G4endl;
	      edep = (*gammaCollection)[ii-1]->GetEdep() - residualEnergy;
	      if( edep <= 0. ) edep = 0.;
	      (*gammaCollection)[ii-1]->SetEdep(edep);
	      //	      G4cout << " Deposited energy now is " << (*gammaCollection)[ii-1]->GetEdep()/keV << G4endl;
	      residualEnergy = 0.;
	    }
	  }
	}
	if( ChiEStato != "gamma" ) return false;
   

	///////////////////////////////////////////////
	///////// Treatment of secondary gammas
	///////////////////////////////////////////////

	padre = aStep->GetTrack()->GetParentID();
	//	G4cout << "ParentID = " << padre << G4endl;

	// LR: Primary gammas from in-beam simulation have ParentID = 2. 
	//     In source simulations, primaries have ParentID = 0.
	if( padre > primaryGammaID ) {
	  ////////////////////////////////////////////////////////////////////////
	  //// Should not count very-low-energy gammas from photoelectric
	  ///  or bremsstrahlung
	  ////////////////////////////////////////////////////////////////////////
	  if( aStep->GetTrack()->GetCreatorProcess() )
	    ilGeneratore = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
	  if( ilGeneratore != "annihil" )
	    return false;

	  ////////////////////////////////////////////////////////////////////////
	  //// Only secondary gammas from annihilation vertexes!!!
	  ////////////////////////////////////////////////////////////////////////
	  if( annihilationHappened ) {
	    G4ThreeVector originOfTrack = aStep->GetTrack()->GetVertexPosition();
	    //	    G4cout << "          originOfTrack is " << originOfTrack/mm << G4endl;
	    if( (positionOfAnnihilation-originOfTrack).mag() > 0.0001*mm ) return false;
	  }
	}
	///////////////////////////////////////////////////////////////////
	/////// Should take care of in-flight annihilations!
	//////////////////////////////////////////////////////////////////
	if( ilProcesso == "conv" || ilProcesso == "LowEnConversion" ) {
	  additionalEnergy = residualEnergy;
	  //	  if( additionalEnergy )
	  //	    G4cout << " additionalEnergy is " << additionalEnergy << G4endl;
	}

	// LR: G4Step->GetDeltaEnergy() has been made obsolete, so now we do it "by hand."
	//	G4cout << std::fixed << std::setprecision(4) << std::setw(10)
	//	       << "aStep->GetDeltaEnergy() = " << aStep->GetDeltaEnergy() 
	//	       << "   deltaEnergy = " << deltaEnergy << G4endl;

	deltaEnergy = aStep->GetPostStepPoint()->GetKineticEnergy() - aStep->GetPreStepPoint()->GetKineticEnergy();
	G4cout << "Post-step kinetic energy = " << aStep->GetPostStepPoint()->GetKineticEnergy()/keV << G4endl;
	if( ilProcesso != "conv" && ilProcesso != "LowEnConversion" )
	  edep  = -1.*deltaEnergy; // This gives the total energy deposited by the gamma and its secondaries.
	else {      // pair production
	  //	  G4cout << "                deltaEnergy = " << deltaEnergy/keV << G4endl;
	  //	  G4cout << " pair production correction = " << (2. * electron_mass_c2 - additionalEnergy)/keV << G4endl;
	  edep  =  -1.*deltaEnergy - 2. * electron_mass_c2 - additionalEnergy;  // E_gamma - 2 m_e
	  if( edep < 0. ) edep = 0.;
	}

	if(edep==0.)
	  return false;

	break;

      case 2:
	// Like case1, but consider all secondary gammas. (This means also considering e+/e- brehmstrahlung.)

	/// Consider only gammas and positrons
	ChiEStato = aStep->GetTrack()->GetDefinition()->GetParticleName();
	ilProcesso = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	//	G4cout << ChiEStato << " interacted through a " << ilProcesso << G4endl;
	/// positrons: just store the position!!!
	if( ChiEStato == "e+" && ilProcesso == "annihil" ) {
	  positionOfAnnihilation = aStep->GetPostStepPoint()->GetPosition();
	  /// some positron decays in flight, should correct for that ...
	  //	  residualEnergy = aStep->GetTrack()->GetKineticEnergy(); // LR: always zero?
	  // LR: This gives us our escape peaks:
	  residualEnergy = -1.*(aStep->GetPostStepPoint()->GetKineticEnergy() 
				- aStep->GetPreStepPoint()->GetKineticEnergy());
	  annihilationHappened = true;
	  //	  G4cout << " e+ annihilated at " << positionOfAnnihilation/mm << G4endl; 
	  //	  G4cout << " e+ residual energy is " << residualEnergy/keV << G4endl;
	  if( residualEnergy ) {
	    //	    G4cout << ChiEStato << " interacted through a " << ilProcesso << G4endl;
	    //	    G4cout << " e+ residual energy is " << residualEnergy/keV << G4endl;
	    G4int ii = gammaCollection->entries();
	    //	    G4cout << " number of entries is " << ii << G4endl;
	    if( ii ) {
	      //	      G4cout << " Deposited energy was " << (*gammaCollection)[ii-1]->GetEdep()/keV << G4endl;
	      edep = (*gammaCollection)[ii-1]->GetEdep() - residualEnergy;
	      if( edep <= 0. ) edep = 0.;
	      (*gammaCollection)[ii-1]->SetEdep(edep);
	      //	      G4cout << " Deposited energy now is " << (*gammaCollection)[ii-1]->GetEdep()/keV << G4endl;
	      residualEnergy = 0.;
	    }
	  }
	}
	if( ChiEStato != "gamma" ) return false;
   

	///////////////////////////////////////////////
	///////// Treatment of secondary gammas
	///////////////////////////////////////////////

	padre = aStep->GetTrack()->GetParentID();
	//	G4cout << "ParentID = " << padre << G4endl;

	// LR: Primary gammas from in-beam simulation have ParentID = 2. 
	//     In source simulations, primaries have ParentID = 0.
	if( padre > primaryGammaID ) {
	  ////////////////////////////////////////////////////////////////////////
	  //// Should not count very-low-energy gammas from photoelectric
	  ///  or bremsstrahlung
	  ////////////////////////////////////////////////////////////////////////
	  if( aStep->GetTrack()->GetCreatorProcess() )
	    ilGeneratore = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
	  if( ilGeneratore != "annihil" )
	    return false;

	  ////////////////////////////////////////////////////////////////////////
	  //// Only secondary gammas from annihilation vertexes!!!
	  ////////////////////////////////////////////////////////////////////////
	  if( annihilationHappened ) {
	    G4ThreeVector originOfTrack = aStep->GetTrack()->GetVertexPosition();
	    //	    G4cout << "          originOfTrack is " << originOfTrack/mm << G4endl;
	    if( (positionOfAnnihilation-originOfTrack).mag() > 0.0001*mm ) return false;
	  }
	}
	///////////////////////////////////////////////////////////////////
	/////// Should take care of in-flight annihilations!
	//////////////////////////////////////////////////////////////////
	if( ilProcesso == "conv" || ilProcesso == "LowEnConversion" ) {
	  additionalEnergy = residualEnergy;
	  //	  if( additionalEnergy )
	  //	    G4cout << " additionalEnergy is " << additionalEnergy << G4endl;
	}

	// LR: G4Step->GetDeltaEnergy() has been made obsolete, so now we do it "by hand."
	//	G4cout << std::fixed << std::setprecision(4) << std::setw(10)
	//	       << "aStep->GetDeltaEnergy() = " << aStep->GetDeltaEnergy() 
	//	       << "   deltaEnergy = " << deltaEnergy << G4endl;

	deltaEnergy = aStep->GetPostStepPoint()->GetKineticEnergy() - aStep->GetPreStepPoint()->GetKineticEnergy();
	G4cout << "Post-step kinetic energy = " << aStep->GetPostStepPoint()->GetKineticEnergy()/keV << G4endl;
	if( ilProcesso != "conv" && ilProcesso != "LowEnConversion" )
	  edep  = -1.*deltaEnergy; // This gives the total energy deposited by the gamma and its secondaries.
	else {      // pair production
	  //	  G4cout << "                deltaEnergy = " << deltaEnergy/keV << G4endl;
	  //	  G4cout << " pair production correction = " << (2. * electron_mass_c2 - additionalEnergy)/keV << G4endl;
	  edep  =  -1.*deltaEnergy - 2. * electron_mass_c2 - additionalEnergy;  // E_gamma - 2 m_e
	  if( edep < 0. ) edep = 0.;
	}

	if(edep==0.)
	  return false;

	break;

      default:

	// Track everything (the primary gamma, all secondary electrons, and all secondary gammas).
	edep = aStep->GetTotalEnergyDeposit();
	if(edep<0.001*eV) return false;

      }

      G4int detCode = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(depth);
      G4int detNum  = detCode%1000;

      G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();

      ////////////////////////////////////////////////////////////////////
      /// Position of the interaction point in the solid reference frame
      ///////////////////////////////////////////////////////////////////
      G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
      G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
      G4VPhysicalVolume* topVolume = theTouchable->GetVolume(depth);
  
      G4ThreeVector frameTrans = topVolume->GetFrameTranslation();
                                  
      const G4RotationMatrix* rot = topVolume->GetFrameRotation();
                                  
      G4RotationMatrix frameRot;
      if( rot )
	frameRot = *rot;
  
      G4ThreeVector posSol = frameRot( position );
      posSol += frameRot( frameTrans );
    
      //////////////////////
      /// Segment number
      //////////////////////
      G4int segCode = 0;

      segCode = theDetector->GetGretina()->GetSegmentNumber( 0, detCode, posSol );

      TrackerGammaHit* newHit = new TrackerGammaHit();

      newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());

      newHit->SetParticleID(aStep->GetTrack()->GetDefinition()->GetParticleName());
      newHit->SetProcess(aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
      newHit->SetParentTrackID(aStep->GetTrack()->GetParentID());
      if( aStep->GetTrack()->GetCreatorProcess() )
	newHit->SetCreatorProcess(aStep->GetTrack()->GetCreatorProcess()->GetProcessName());
      else
	newHit->SetCreatorProcess("primary");
      newHit->SetDetNumb    (detNum);
      newHit->SetSegNumb    (segCode);
      newHit->SetEdep       (edep);
      newHit->SetPos        (position);
      newHit->SetTrackOrigin(aStep->GetTrack()->GetVertexPosition());

      gammaCollection->insert( newHit );
      newHit->Draw();
      //  getc(stdin);
  //   newHit->Print(); 
      return true;
      }
  else
    {
      //    G4cout << "Energy deposit in the " << name << G4endl;
      //    G4cout << "E="<<G4BestUnit(edep,"Energy")<<G4endl;
      //    G4cout <<"Event ignored" <<G4endl;
	    //	    getc(stdin);
      return false;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackerGammaSD::EndOfEvent(G4HCofThisEvent* HCE)
{
   G4int i;
   G4int NbHits = gammaCollection->entries();

        if (NbHits>0&&print) 
	  { 

	    G4RunManager* runManager = G4RunManager::GetRunManager();
	    
	    G4cout << "\n--------> event " << runManager->GetCurrentEvent()->GetEventID() << ", "
		   << NbHits << " hits for gamma tracking: " << G4endl;

	    for (i=0;i<NbHits;i++) (*gammaCollection)[i]->Print();

	  }


  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, gammaCollection ); 
 }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

