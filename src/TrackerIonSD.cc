#include "TrackerIonSD.hh"

//--------------------------------------------------------------------
TrackerIonSD::TrackerIonSD(G4String name)
  :G4VSensitiveDetector(name)
{
  G4String HCname;
  collectionName.insert(HCname="ionCollection");
  print=false; //LR (formerly not initialized)
 
}

//--------------------------------------------------------------------

  TrackerIonSD::~TrackerIonSD()
{ 
  
}

//--------------------------------------------------------------------

void TrackerIonSD::Initialize(G4HCofThisEvent*)
{ 
 

    ionCollection = new TrackerIonHitsCollection
                          (SensitiveDetectorName,collectionName[0]);  
  
}
//--------------------------------------------------------------------


G4bool TrackerIonSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{

 G4StepPoint*   vi;
 G4StepPoint*   vf;

 const G4DynamicParticle* aParticle= aStep->GetTrack()->GetDynamicParticle();
 const G4String type =  aParticle->GetDefinition()->GetParticleType();
 const G4int    trackID = aStep->GetTrack()->GetTrackID(); 
 const G4String name =  aParticle->GetDefinition()->GetParticleName();
 const G4double len=aStep->GetStepLength();
 // G4cout<<G4endl;

  if(type=="nucleus")
   {     

     vi=aStep->GetPreStepPoint();    
     G4String vname=vi->GetPhysicalVolume()->GetName();
     if(len>0.)
     {
         TrackerIonHit* newIonHitI= new TrackerIonHit();           
	 newIonHitI->SetTrackID(trackID);
	 newIonHitI->SetParticleID(name);
	 newIonHitI->SetBeta(vi->GetBeta());
	 newIonHitI->SetTheta(vi->GetMomentumDirection().theta());
	 newIonHitI->SetPhi(vi->GetMomentumDirection().phi());
	 newIonHitI->SetTime(vi->GetProperTime());
	 newIonHitI->SetLabTime(vi->GetLocalTime());
	 newIonHitI->SetGlobTime(vi->GetGlobalTime());
	 newIonHitI->SetKE(vi->GetKineticEnergy());
	 newIonHitI->SetPos(vi->GetPosition());
	 newIonHitI->SetVolName(vname);
	 //	 G4cout<<" Edep I "<<aStep->GetTotalEnergyDeposit()<<G4endl;
	 newIonHitI->SetEdep(aStep->GetTotalEnergyDeposit());
	 newIonHitI->SetLength(len);
	 newIonHitI->Draw();
	 ionCollection->insert(newIonHitI);
       }
     G4TrackStatus TrackStatus;
     TrackStatus=aStep->GetTrack()->GetTrackStatus();

     if(TrackStatus==fStopButAlive||TrackStatus==fStopAndKill)
       {
	 TrackerIonHit* newIonHitF= new TrackerIonHit();   	
	 vf=aStep->GetPostStepPoint();

	 if(vf->GetProcessDefinedStep()->GetProcessName()=="Decay")
	   newIonHitF->SetDecayFlag();

	 if(vf->GetProcessDefinedStep()->GetProcessName()=="Reaction")
	   newIonHitF->SetReactionFlag();

	 newIonHitF->SetTrackID(trackID);
	 newIonHitF->SetParticleID(name);
	 newIonHitF->SetBeta(vf->GetBeta());
	 newIonHitF->SetTheta(vf->GetMomentumDirection().theta());
	 newIonHitF->SetPhi(vf->GetMomentumDirection().phi());
	 newIonHitF->SetTime(vf->GetProperTime());
	 newIonHitF->SetLabTime(vf->GetLocalTime());
	 newIonHitF->SetGlobTime(vf->GetGlobalTime());
	 newIonHitF->SetKE(vf->GetKineticEnergy());
	 newIonHitF->SetPos(vf->GetPosition());
	 newIonHitF->SetVolName(vname);	
	 newIonHitF->SetEdep(aStep->GetTotalEnergyDeposit());
	 //	 G4cout<<" Edep F "<<aStep->GetTotalEnergyDeposit()<<G4endl;
	 newIonHitF->SetLength(0.);
	 newIonHitF->Draw();
	 ionCollection->insert(newIonHitF);

       }

   
   }

 
  
  return true;
}
//--------------------------------------------------------------------

void TrackerIonSD::EndOfEvent(G4HCofThisEvent* HCE)
{

  G4int i;
  G4int NbHits = ionCollection->entries();  

  if(NbHits>0)
    {
 
      //       (*ionCollection)[0]->SetGunFlag();

       // for (i=0;i<NbHits-1;i++)
       // 	 if ((*ionCollection)[i]->GetVolName()=="expHall"&&(*ionCollection)[i+1]->GetVolName()=="target")
       // 	     (*ionCollection)[i+1]->SetTargetFaceFlag();


       // for (i=1;i<NbHits;i++) 
       // 	 if ((*ionCollection)[i-1]->GetVolName()=="target"&&(*ionCollection)[i]->GetVolName()=="expHall")
       // 	     (*ionCollection)[i]->SetTargetBackFlag();

     if (print) 
     {	
       std::ios init(NULL);
       init.copyfmt(G4cout);

       G4cout<<G4endl;
       G4cout << "-------->Hits Collection: in this event there are "
	      << NbHits
	      << " hits for ion tracking: " << G4endl;

       G4cout << std::setw(2) << " " << "ID" <<" "
	      << std::setw(15)<< "PID"<<" "
	      << std::setprecision(4)<<std::setw(8)<<std::fixed
	      << "KE [GeV]"<<" "<<std::setprecision(6)<<std::setw(8)
	      << "beta"<< " "
	      <<std::setprecision(4)<<std::setw(10)<<std::right
	      <<"th[mrad]"<<" "<<std::setw(12)<<std::right
	      <<"phi [mrad]"<<" "<<std::setw(10)<<std::right	 
	      <<"X [mm]" <<" "<<std::setw(10)<<std::right
	      <<"Y [mm]"<<" "<<std::setw(10)<<std::right
	      <<"Z [mm]"<<" "<<std::setw(12)<<std::right
	      <<"tau [ps]"<<" "<<std::setw(12)<<std::right
	      <<"t [ps]"<<" "<<std::setw(12)<<std::right
	      <<"T [ps]"<<" "<<std::setw(10)<<std::right
	      <<"Ed [MeV]"<<" " <<std::setw(10)<<std::right
	      <<"Length [mm]"
	      << G4endl;
       
       for (i=0;i<NbHits;i++) (*ionCollection)[i]->Print();
	   
       G4cout.copyfmt( init );
     
     }
     
    }  
 
 static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, ionCollection ); 
 }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

