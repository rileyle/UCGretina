#include "PhysicsList.hh"

PhysicsList::PhysicsList()
{
  BeamOut = NULL;
}

PhysicsList::PhysicsList(Outgoing_Beam* BO):BeamOut(BO)
{;}

PhysicsList::~PhysicsList()
{;}

void PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  //  ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
  
}

void PhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
  ConstructEM();
}

void PhysicsList::ConstructEM()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    G4cout<<"++++ Particle name ="<<particleName<<G4endl;
    if (particleName == "gamma") {

      // Livermore Low Energy EM processes ====================================
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePhotoElectricModel();
      thePhotoElectricEffect->SetEmModel(theLivermorePhotoElectricModel);
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermoreComptonModel* theLivermoreComptonModel = new G4LivermoreComptonModel();
      theComptonScattering->SetEmModel(theLivermoreComptonModel);
      pmanager->AddDiscreteProcess(theComptonScattering);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermoreGammaConversionModel();
      theGammaConversion->SetEmModel(theLivermoreGammaConversionModel);
      pmanager->AddDiscreteProcess(theGammaConversion);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
      theRayleigh->SetEmModel(theRayleighModel);
      pmanager->AddDiscreteProcess(theRayleigh);

    }
    else if (particleName == "e-") {

      pmanager->AddProcess(new G4eMultipleScattering,-1,1,1);

      // Livermore Low Energy EM processes =====================================
      G4eIonisation* eIoni = new G4eIonisation();
      G4LivermoreIonisationModel* theLivermoreIonisationModel = new G4LivermoreIonisationModel();
      eIoni->SetEmModel(theLivermoreIonisationModel);
      eIoni->SetStepFunction(0.2, 100*um); //     
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      eBrem->SetEmModel(new G4LivermoreBremsstrahlungModel());
      pmanager->AddProcess(eBrem,         -1,-3, 3);

    }
    else if (particleName == "e+") {

      pmanager->AddProcess(new G4eMultipleScattering,-1,1,1);

      // Livermore Low Energy EM models: no e+ (used standard in testing). =====
      // Standard EM processes =================================================
      pmanager->AddProcess(new G4eIonisation,       -1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,3,3);     
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4); 

    }
    else if( particleName == "GenericIon" ) {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2, 2);
      pmanager->AddProcess(new Reaction(BeamOut),    -1,-1, 3);
      pmanager->AddProcess(new G4StepLimiter,        -1,-1, 4);

    } 

  }

}

void PhysicsList::SetCuts()
{
  // uppress error messages even in case e/gamma/proton do not exist            
  G4int temp = GetVerboseLevel();                                                
  SetVerboseLevel(0);                                                           
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   

  // Retrieve verbose level
  SetVerboseLevel(temp);  
}
