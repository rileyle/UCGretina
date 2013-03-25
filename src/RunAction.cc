#include "RunAction.hh"

RunAction::RunAction(DetectorConstruction* detector, Outgoing_Beam* BO,EventAction* ev): myDetector(detector), BeamOut(BO),evaction(ev)
{
  
}


RunAction::~RunAction()
{

}

void RunAction::BeginOfRunAction(const G4Run*)
{

  G4cout<<" Begining of run "<<G4endl;
  evaction->openEvfile(); //TB
  if(BeamOut->ReactionOn()) 
    G4cout<<" Simulating scattered ions "<<G4endl;
  else
     G4cout<<" Simulating unscattered ions "<<G4endl;
//  myDetector->GetSeGA()->CalcSegCenters(BeamOut->GetPosDopp());       // LR
//   results->resetTauAverage();
//   results->resetKEAverage();

}


 
void RunAction::EndOfRunAction(const G4Run*)
{
	evaction->closeEvfile();
  G4cout<<G4endl;
//  G4cout<<" End of run "<<G4endl;
 // G4cout<<" Average decay time:"<<G4endl;
//  G4cout<<"                set: "<<G4BestUnit(BeamOut->getTime(),"Time")<<G4endl;
//  G4cout<<"          simulated: "<<G4BestUnit(results->getTauAverage(),"Time")<<G4endl;

//  G4cout<<" Average energy deposit in the target"<<G4endl;
//  G4cout<<"          simulated: "<<std::setw(9)<<std::setprecision(3)<<results->getTarDepAverage()/MeV<<" +/- "<<std::setw(4)<<std::setprecision(2)<<results->getSigmaTarDepAverage()/MeV<<" MeV "<<G4endl;



//  G4cout<<" Average Kinetic Energy behind the target"<<G4endl;
//  G4cout<<"          simulated: "<<std::setw(9)<<std::setprecision(7)<<results->getKEAverage()/GeV<<" +/- "<<std::setw(4)<<std::setprecision(3)<<results->getSigmaKEAverage()/GeV<<" GeV "<<G4endl;
 
}

