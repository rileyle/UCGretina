//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm1/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
// 
// $Id: PhysicsList.cc 100290 2016-10-17 08:47:55Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsList_Messenger.hh"
 
//#include "PhysListEmStandard.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"

#include "G4PolarizedPhotoElectricEffect.hh"
#include "G4PolarizedCompton.hh"
#include "G4PolarizedGammaConversion.hh"

#include "DetectorConstruction.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Gamma.hh"

#ifdef NEUTRONS
// From LBE for neutrons
// Builder for all stopping processes
#include "G4StoppingPhysics.hh"
#include "G4MaxTimeCuts.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(DetectorConstruction* det) 
  : G4VModularPhysicsList(), fDet(det)
{
  theMessenger = new PhysicsList_Messenger(this);
  SetVerboseLevel(1);

  // EM physics
  //  fEmName = G4String("emstandard_opt4_mod");
  fEmName = G4String("emstandard_opt4");
  fEmPhysicsList = new G4EmStandardPhysics_option4();
  
  G4LossTableManager::Instance();
  // fix lower limit for cut
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(10*eV, 1*GeV);
  SetDefaultCutValue(1*mm);

  BeamOut = NULL;

  usePolar = false;
  
#ifdef NEUTRONS
  // From LBE for neutrons
  stoppingPhysics = new G4StoppingPhysics;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{
  delete theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
    G4BosonConstructor  pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  // Transportation
  //
  AddTransportation();

#ifdef NEUTRONS  
  // From LBE for neutrons
  auto myParticleIterator=G4ParticleTable::GetParticleTable()->GetIterator();
  myParticleIterator->reset();
  while( (*(myParticleIterator))() ){
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    // time cuts for ONLY neutrons:
    if(particleName == "neutron") 
      pmanager->AddDiscreteProcess(new G4MaxTimeCuts());
    // Energy cuts to kill charged (embedded in method) particles:
    //    pmanager->AddDiscreteProcess(new G4MinEkineCuts()); // LR: need this?
  }	

  // From LBE for neutrons
  ConstructHad();
#endif
  
  // Electromagnetic physics list
  //
  fEmPhysicsList->ConstructProcess();
  
  if(usePolar){
    G4ProcessManager *gpMan = G4Gamma::Gamma()->GetProcessManager();
    G4ProcessVector* pv = gpMan->GetProcessList();
    for(unsigned int i=0;i<pv->entries();i++){
      if((*pv)[i]->GetProcessName()=="phot"){
	gpMan->RemoveProcess((*pv)[i]);
	gpMan->AddDiscreteProcess(new G4PolarizedPhotoElectricEffect);
      }
      if((*pv)[i]->GetProcessName()=="compt"){
	gpMan->RemoveProcess((*pv)[i]);
	gpMan->AddDiscreteProcess(new G4PolarizedCompton());
      }
      if((*pv)[i]->GetProcessName()=="conv"){
	gpMan->RemoveProcess((*pv)[i]);
	gpMan->AddDiscreteProcess(new G4PolarizedGammaConversion);
      }
    }
  }

  // Beam Reaction
  AddReaction();

  // Decay Process
  //
  //  AddDecay();
    
  // Decay Process
  //
  AddRadioactiveDecay();  

  // step limitation (as a full process)
  //  
  //  AddStepMax();    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>1) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  
  if (name == fEmName) return;

  //if (name == "local") {
  //  fEmName = name;
  //  delete fEmPhysicsList;
  //  fEmPhysicsList = new PhysListEmStandard(name);
    
  // } else 

  if (name == "emstandard_opt0") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics();

  } else if (name == "emstandard_opt1") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option1();

  } else if (name == "emstandard_opt2") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "emstandard_opt3") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option3();
    
  } else if (name == "emstandard_opt4") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysics_option4();
        
  } else if (name == "emstandardSS") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsSS();
    
  } else if (name == "emstandardGS") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsGS();
    
  } else if (name == "emstandardWVI") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmStandardPhysicsWVI();
    
  } else if (name == "emlivermore") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLivermorePhysics();

  } else if (name == "empenelope") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmPenelopePhysics();
            
  } else if (name == "emlowenergy") {
    fEmName = name;
    delete fEmPhysicsList;
    fEmPhysicsList = new G4EmLowEPPhysics();
            
  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhysicsListHelper.hh"
#include "Reaction.hh"
#include "G4StepLimiter.hh"

void PhysicsList::AddReaction()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    
  Reaction* react = new Reaction(BeamOut);

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();
    G4String particleType = particle->GetParticleType();
    if ( particleType == "nucleus" )
      ph->RegisterProcess(react, particle);    
    ph->RegisterProcess(new G4StepLimiter, particle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhysicsListHelper.hh"
#include "G4Decay.hh"

void PhysicsList::AddDecay()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    
  // Decay Process
  //
  G4Decay* fDecayProcess = new G4Decay();

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    if (fDecayProcess->IsApplicable(*particle)) 
      ph->RegisterProcess(fDecayProcess, particle);    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"
#include "G4NuclideTable.hh"

void PhysicsList::AddRadioactiveDecay()
{  
  G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
  
  radioactiveDecay->SetARM(true);                //Atomic Rearangement
  
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();  
  ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

  // mandatory for G4NuclideTable
  //
  G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// #include "G4ProcessManager.hh"
// #include "StepMax.hh"

// void PhysicsList::AddStepMax()
// {
//   // Step limitation seen as a process
//   StepMax* stepMaxProcess = new StepMax();

//   auto particleIterator=GetParticleIterator();
//   particleIterator->reset();
//   while ((*particleIterator)()){
//       G4ParticleDefinition* particle = particleIterator->value();
//       G4ProcessManager* pmanager = particle->GetProcessManager();

//       if (stepMaxProcess->IsApplicable(*particle))
//         {
//           pmanager ->AddDiscreteProcess(stepMaxProcess);
//         }
//   }
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Material.hh"

void PhysicsList::GetRange(G4double val)
{
  //  G4LogicalVolume* lBox = fDet->GetWorld()->GetLogicalVolume();
  G4LogicalVolume* lBox = fDet->HallLog();
  G4ParticleTable* particleTable =  G4ParticleTable::GetParticleTable();
  const G4MaterialCutsCouple* couple = lBox->GetMaterialCutsCouple();
  const G4Material* currMat = lBox->GetMaterial();

  G4ParticleDefinition* part;
  G4double cut;
  part = particleTable->FindParticle("e-");
  cut = G4LossTableManager::Instance()->GetRange(part,val,couple);
  G4cout << "material : " << currMat->GetName()       << G4endl;
  G4cout << "particle : " << part->GetParticleName()  << G4endl;
  G4cout << "energy   : " << G4BestUnit(val,"Energy") << G4endl;
  G4cout << "range    : " << G4BestUnit(cut,"Length") << G4endl;
}

void PhysicsList::SetGammaAngularCorrelations(bool val){
  G4NuclearLevelData::GetInstance()->GetParameters()->SetCorrelatedGamma(val);
  G4cout<<"Set correlated gamma "<<val<<G4endl;
}

void PhysicsList::SetUsePolarizedPhysics(bool use){
  usePolar=use;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef NEUTRONS
// LR: From LBE for neutrons

// Hadronic processes ////////////////////////////////////////////////////////

#include "G4HadronicParameters.hh"

// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4AntiNuclElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsNeutronElasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"  
#include "G4ChipsKaonMinusElasticXS.hh"
#include "G4ChipsKaonPlusElasticXS.hh"
#include "G4ChipsKaonZeroElasticXS.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4CrossSectionElastic.hh"

// Inelastic processes:
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// FTFP + BERT model
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGNuclNuclXsc.hh"

// Neutron high-precision models: <20 MeV
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronRadCapture.hh"

// Binary light ion cascade for alpha, deuteron and triton
#include "G4BinaryLightIonReaction.hh"

// ConstructHad()
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering and Inelastic scattering.
// F.W.Jones  09-JUL-1998
 void PhysicsList::ConstructHad() 
{
  // Elastic scattering
  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE(); 

  const G4double elastic_elimitAntiNuc = 100.0*CLHEP::MeV;
  G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
  elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
  G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
  G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
  elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );

  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*CLHEP::GeV;
  const G4double theFTFMin1 =    4.0*CLHEP::GeV;
  const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  const G4double theBERTMin0 =   0.0*CLHEP::GeV;
  const G4double theBERTMin1 =  19.0*CLHEP::MeV;
  const G4double theBERTMax =    5.0*CLHEP::GeV;
  const G4double theHPMin =      0.0*CLHEP::GeV;
  const G4double theHPMax =     20.0*CLHEP::MeV;
  const G4double theIonBCMin =   0.0*CLHEP::GeV;
  const G4double theIonBCMax =   5.0*CLHEP::GeV;


  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  // Binary Cascade
  G4BinaryLightIonReaction * theIonBC = new G4BinaryLightIonReaction( thePreEquilib );
  theIonBC->SetMinEnergy( theIonBCMin );
  theIonBC->SetMaxEnergy( theIonBCMax );

  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4ComponentGGNuclNuclXsc * ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet * theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);

  auto myParticleIterator=G4ParticleTable::GetParticleTable()->GetIterator();
  myParticleIterator->reset();
  while ((*(myParticleIterator))()) 
    {
      G4ParticleDefinition* particle = myParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4PionPlusInelasticProcess* theInelasticProcess = new G4PionPlusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( new G4BGGPionInelasticXS( G4PionPlus::Definition() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	} 

      else if (particleName == "pi-") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4PionMinusInelasticProcess* theInelasticProcess = new G4PionMinusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( new G4BGGPionInelasticXS( G4PionMinus::Definition() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "kaon+") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusElasticXS::Default_Name()));
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonPlusInelasticProcess* theInelasticProcess = new G4KaonPlusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "kaon0S") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroElasticXS::Default_Name()));
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonZeroSInelasticProcess* theInelasticProcess = new G4KaonZeroSInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "kaon0L") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroElasticXS::Default_Name()));
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonZeroLInelasticProcess* theInelasticProcess = new G4KaonZeroLInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 ); 
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "kaon-") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusElasticXS::Default_Name()));
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonMinusInelasticProcess* theInelasticProcess = new G4KaonMinusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "proton") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->AddDataSet( new G4BGGNucleonElasticXS( G4Proton::Proton() ) );
          theElasticProcess->RegisterMe( elastic_chip );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4ProtonInelasticProcess* theInelasticProcess =  new G4ProtonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "anti_proton") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4AntiProtonInelasticProcess* theInelasticProcess = new G4AntiProtonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "neutron") {
	// elastic scattering
	G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
	elastic_neutronChipsModel->SetMinEnergy( 19.0*CLHEP::MeV );
        theElasticProcess->RegisterMe( elastic_neutronChipsModel );
	G4ParticleHPElastic * theElasticNeutronHP = new G4ParticleHPElastic;
        theElasticNeutronHP->SetMinEnergy( theHPMin );
        theElasticNeutronHP->SetMaxEnergy( theHPMax );
	theElasticProcess->RegisterMe( theElasticNeutronHP );
	theElasticProcess->AddDataSet( new G4ParticleHPElasticData );
	pmanager->AddDiscreteProcess( theElasticProcess );
	// inelastic scattering
	G4NeutronInelasticProcess* theInelasticProcess = new G4NeutronInelasticProcess("inelastic");
        theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
	theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel1 );
	G4ParticleHPInelastic * theNeutronInelasticHPModel = new G4ParticleHPInelastic;
        theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
        theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
	theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
	theInelasticProcess->AddDataSet( new G4ParticleHPInelasticData );
	pmanager->AddDiscreteProcess(theInelasticProcess);
	// capture
	G4HadronCaptureProcess* theCaptureProcess = new G4HadronCaptureProcess;
	G4ParticleHPCapture * theNeutronCaptureHPModel = new G4ParticleHPCapture;
        theNeutronCaptureHPModel->SetMinEnergy( theHPMin );
        theNeutronCaptureHPModel->SetMaxEnergy( theHPMax );
	G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
 	theNeutronRadCapture->SetMinEnergy(theHPMax*0.99); 
	theCaptureProcess->RegisterMe( theNeutronCaptureHPModel );
	theCaptureProcess->RegisterMe( theNeutronRadCapture);
	theCaptureProcess->AddDataSet( new G4ParticleHPCaptureData );
	theCaptureProcess->AddDataSet((G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name()));
	pmanager->AddDiscreteProcess(theCaptureProcess);
      }
      else if (particleName == "anti_neutron") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4AntiNeutronInelasticProcess* theInelasticProcess = new G4AntiNeutronInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "deuteron") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( theGGNuclNuclData );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4DeuteronInelasticProcess* theInelasticProcess = new G4DeuteronInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
	  theInelasticProcess->RegisterMe( theIonBC );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "triton") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( theGGNuclNuclData );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4TritonInelasticProcess* theInelasticProcess = new G4TritonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
	  theInelasticProcess->RegisterMe( theIonBC );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "alpha") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( theGGNuclNuclData );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4AlphaInelasticProcess* theInelasticProcess = new G4AlphaInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
	  theInelasticProcess->RegisterMe( theIonBC );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
    }	// while ((*(myParticleIterator))()) 

  // Add stopping processes with builder
  stoppingPhysics->ConstructProcess();
}
#endif
