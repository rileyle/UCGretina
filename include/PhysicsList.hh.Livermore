#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"

//gamma
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"

//e-/e+
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4eplusAnnihilation.hh"

// Livermore Low Energy EM Models
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// Penelope Low Energy EM Models
//#include "G4PenelopeIonisationModel.hh"
//#include "G4PenelopeBremsstrahlungModel.hh"
//#include "G4PenelopeAnnihilationModel.hh"
//#include "G4PenelopeComptonModel.hh"
//#include "G4PenelopePhotoElectricModel.hh"
//#include "G4PenelopeGammaConversionModel.hh"
//#include "G4PenelopeRayleighModel.hh"

// ions
#include "G4ionIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4IonTable.hh"
#include "G4IonConstructor.hh"
#include "Reaction.hh"
#include "Outgoing_Beam.hh"
#include "G4StepLimiter.hh"

class PhysicsList: public G4VUserPhysicsList
{
  public:
  PhysicsList();
    PhysicsList(Outgoing_Beam*);
    ~PhysicsList();

  void SetOutgoingBeam(Outgoing_Beam *BO) {BeamOut = BO;}
  protected:
    // Construct particle and physics process
  void ConstructParticle();
  void ConstructProcess();
  void ConstructEM();
  void SetCuts();

private:
  Outgoing_Beam* BeamOut;
  G4double  cutForGamma;      //LR
  G4double  cutForElectron;   //LR
  G4double  cutForPositron;   //LR

};

#endif
