#ifdef POL
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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4ITDecay.cc                                                      //
//  Author: D.H. Wright (SLAC)                                                //
//  Date:   14 November 2014                                                  //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "G4ITDecay.hh"
#include "G4IonTable.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4PhotonEvaporation.hh"
#include "G4RadioactiveDecay.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShells.hh"
#include "G4Electron.hh"
#include "G4LossTableManager.hh"
#include "G4Fragment.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PolarizationTransition.hh"
#include "G4Clebsch.hh"
#include "G4NuclearPolarizationStore.hh"


G4ITDecay::G4ITDecay(const G4ParticleDefinition* theParentNucleus,
                     const G4double& branch, const G4double& Qvalue,
                     const G4double& excitationE, G4PhotonEvaporation* aPhotoEvap)
 : G4NuclearDecay("IT decay", IT, excitationE, noFloat), transitionQ(Qvalue), 
   applyARM(true), photonEvaporation(aPhotoEvap)
{
  SetParent(theParentNucleus);  // Store name of parent nucleus, delete G4MT_parent 
  SetBR(branch);

  parentZ = theParentNucleus->GetAtomicNumber();
  parentA = theParentNucleus->GetAtomicMass(); 

  SetNumberOfDaughters(1);
  G4IonTable* theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  SetDaughter(0, theIonTable->GetIon(parentZ, parentA, excitationE, noFloat) );
}


G4ITDecay::~G4ITDecay()
{}


G4DecayProducts* G4ITDecay::DecayIt(G4double)
{
  // Fill G4MT_parent with theParentNucleus (stored by SetParent in ctor)  
  CheckAndFillParent();

  // Set up final state
  // parentParticle is set at rest here because boost with correct momentum 
  // is done later
  G4LorentzVector atRest(G4MT_parent->GetPDGMass(),
                         G4ThreeVector(0.,0.,0.) );
  G4DynamicParticle parentParticle(G4MT_parent, atRest);
  G4DecayProducts* products = new G4DecayProducts(parentParticle);

  // Let G4PhotonEvaporation do the decay
  G4Fragment parentNucleus(parentA, parentZ, atRest);
  G4double predecayEnergy = parentNucleus.GetExcitationEnergy();

  //Get the polarization of the initial state before it is updated by decay
  std::vector<std::vector<G4complex>> parentPol
    = G4NuclearPolarizationStore::GetInstance()->FindOrBuild(parentNucleus.GetZ_asInt(),parentNucleus.GetA_asInt(),parentNucleus.GetExcitationEnergy())->GetPolarization();
  //parentNucleus.GetNuclearPolarization()->GetPolarization();

  G4Fragment* eOrGamma = photonEvaporation->EmittedFragment(&parentNucleus);
  // Modified nuclide is returned as dynDaughter
  G4IonTable* theIonTable =
    (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable() );
  G4ParticleDefinition* daughterIon =
    theIonTable->GetIon(parentZ, parentA, parentNucleus.GetExcitationEnergy(), 
                        G4Ions::FloatLevelBase(parentNucleus.GetFloatingLevelNumber()));
  G4DynamicParticle* dynDaughter = new G4DynamicParticle(daughterIon,
                                                         parentNucleus.GetMomentum());

  if (eOrGamma) {
    G4DynamicParticle* eOrGammaDyn =
      new G4DynamicParticle(eOrGamma->GetParticleDefinition(),
                            eOrGamma->GetMomentum() );
    eOrGammaDyn->SetProperTime(eOrGamma->GetCreationTime() );
    //If daughter is gamma, calculate its Stokes vector and set it.
    //We don't handle the case of electrons right now.
    if(eOrGamma->GetParticleDefinition()->GetParticleName()==G4String("gamma")){
      G4ThreeVector gammaDir = eOrGammaDyn->GetMomentumDirection();
      //gammaDir.setTheta(CLHEP::pi/2.);
      const G4LevelManager *lman = G4NuclearLevelData::GetInstance()->GetLevelManager(parentZ,parentA);
      size_t pIndex = lman->NearestLevelIndex(predecayEnergy);
      //      size_t dIndex = lman->NearestLevelIndex(parentNucleus.GetExcitationEnergy()-eOrGamma->GetMomentum().getT());
      size_t dIndex = lman->NearestLevelIndex(predecayEnergy-eOrGamma->GetMomentum().getT());
      const G4NucLevel *level = lman->GetLevel(pIndex);
      // Look up the index of the transition.
      size_t tIndex;
      for(tIndex = 0; tIndex< lman->NumberOfTransitions(); tIndex++)
	if(level->FinalExcitationIndex(tIndex) == dIndex) break;
      //      G4double mpRatio = level->MultipolarityRatio(dIndex);
      G4double mpRatio = level->MultipolarityRatio(tIndex);
      G4int JP1=lman->SpinTwo(pIndex);
      G4int JP2=lman->SpinTwo(dIndex);
      //      G4int MP = level->TransitionType(dIndex);
      G4int MP = level->TransitionType(tIndex);
      int Lbar;
      if(MP<99)//Not mixed transition
	Lbar=MP/2;
      else
	Lbar=MP/200;
      G4double P1 = 0.;
      G4double P2 = 0.;//Stokes parameters
      G4double P3 = 0.;
      G4double norm = 0.;//Normalization of Stokes parameters
      G4PolarizationTransition polT;//Just for utilities
      for(size_t k=0;k<parentPol.size();k++){
	G4double pre1=0.;//Precompute the sum over kappa
	G4double pre2=0.;
	G4double pre3=0.;
	G4double preNorm=0.;
	if(k%2){//P3 is the only term which depends on odd-k
	  for(size_t kappa=0;kappa<parentPol[k].size();kappa++){
	    pre3 += (kappa==0?1.:2.)
	      *(G4complex(cos(kappa*gammaDir.phi()),sin(kappa*gammaDir.phi()))
		*G4Clebsch::WignerLittleD(2*k,0,2*kappa,gammaDir.cosTheta())
		*parentPol[k][kappa]).real();
	  }
	  if(MP<99)
	    P3 += 2*sqrt(2*k+1)*pre3*polT.FCoefficient(k,Lbar,Lbar,JP2,JP1);
	  else{
	    P3 += 2*sqrt(2*k+1)*pre3*
	    (polT.FCoefficient(k,Lbar,Lbar,JP2,JP1)
	     +2*mpRatio*polT.FCoefficient(k,Lbar,Lbar+1,JP2,JP1)
	     +mpRatio*mpRatio*polT.FCoefficient(k,Lbar+1,Lbar+1,JP2,JP1)
	     );
	  }
	}
	else{
	  for(size_t kappa=0;kappa<parentPol[k].size();kappa++){
	    preNorm += (kappa==0?1.:2.)
	      *(
		parentPol[k][kappa]
		*G4complex(cos(kappa*gammaDir.phi()),sin(kappa*gammaDir.phi()))
		*G4Clebsch::WignerLittleD(2*k,0,2*kappa,gammaDir.cosTheta())
		).real();
	    pre1 += (kappa==0?1.:2.)
	      *(
		parentPol[k][kappa]
		*(
		  G4complex(cos(kappa*gammaDir.phi()),sin(kappa*gammaDir.phi()))
		  *G4Clebsch::WignerLittleD(2*k,4,2*kappa,gammaDir.cosTheta())
		  //*G4complex(cos(CLHEP::pi),sin(CLHEP::pi))
		  +
		  G4complex(cos(kappa*gammaDir.phi()),sin(kappa*gammaDir.phi()))
		  *G4Clebsch::WignerLittleD(2*k,-4,2*kappa,gammaDir.cosTheta())
		  //*G4complex(cos(CLHEP::pi),-sin(CLHEP::pi))
		  )
		).real();
	    pre2 += (kappa==0?1.:2.)
	      *(
		parentPol[k][kappa]
		*(
		  G4complex(cos(kappa*gammaDir.phi()),sin(kappa*gammaDir.phi()))
		  *G4Clebsch::WignerLittleD(2*k,4,2*kappa,gammaDir.cosTheta())
		  //*G4complex(cos(CLHEP::pi/2.),sin(CLHEP::pi/2.))
		  -
		  G4complex(cos(kappa*gammaDir.phi()),sin(kappa*gammaDir.phi()))
		  *G4Clebsch::WignerLittleD(2*k,-4,2*kappa,gammaDir.cosTheta())
		  //*G4complex(cos(CLHEP::pi/2.),-sin(CLHEP::pi/2.))
		  )
		).imag();
	  }
	  if(MP<99){
	    norm += 2*sqrt(2*k+1)*preNorm*polT.FCoefficient(k,Lbar,Lbar,JP2,JP1);
	    P1 += sqrt(2*k+1)*pre1
	      *(
		(
		 G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,2,2*k)
		 /G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,-2,2*k)
		 *polT.FCoefficient(k,Lbar,Lbar,JP2,JP1)
		 )
	       );
	    P2 += sqrt(2*k+1)*pre2
	      *(
		(
		 G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,2,2*k)
		 /G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,-2,2*k)
		 *polT.FCoefficient(k,Lbar,Lbar,JP2,JP1)
		 )
	       );
	  }
	  else{
	    norm += 2*sqrt(2*k+1)*preNorm*
	    (polT.FCoefficient(k,Lbar,Lbar,JP2,JP1)
	     +2*mpRatio*polT.FCoefficient(k,Lbar,Lbar+1,JP2,JP1)
	     +mpRatio*mpRatio*polT.FCoefficient(k,Lbar+1,Lbar+1,JP2,JP1)
	     );
	    P1 += sqrt(2*k+1)*pre1
	      *(
		(
		 G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,2,2*k)
		 /G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,-2,2*k)
		 *polT.FCoefficient(k,Lbar,Lbar,JP2,JP1)
		 )
		-(
		 2*mpRatio*
		 G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*(Lbar+1),2,2*k)
		 /G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*(Lbar+1),-2,2*k)
		 *polT.FCoefficient(k,Lbar,Lbar+1,JP2,JP1)
		 )
		-(
		  mpRatio*mpRatio*
		  G4Clebsch::ClebschGordanCoeff(2*(Lbar+1),2,2*(Lbar+1),2,2*k)
		  /G4Clebsch::ClebschGordanCoeff(2*(Lbar+1),2,2*(Lbar+1),-2,2*k)
		  *polT.FCoefficient(k,Lbar+1,Lbar+1,JP2,JP1)
		  )
		);
	    P2 += sqrt(2*k+1)*pre2
	      *(
		(
		 G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,2,2*k)
		 /G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,-2,2*k)
		 *polT.FCoefficient(k,Lbar,Lbar,JP2,JP1)
		 )
		-(
		  2*mpRatio*
		  G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*(Lbar+1),2,2*k)
		  /G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*(Lbar+1),-2,2*k)
		  *polT.FCoefficient(k,Lbar,Lbar+1,JP2,JP1)
		  )
		-(
		  mpRatio*mpRatio*
		  G4Clebsch::ClebschGordanCoeff(2*(Lbar+1),2,2*(Lbar+1),2,2*k)
		  /G4Clebsch::ClebschGordanCoeff(2*(Lbar+1),2,2*(Lbar+1),-2,2*k)
		  *polT.FCoefficient(k,Lbar+1,Lbar+1,JP2,JP1)
		  )
		);
	  }
	}
      }
      P1/=(((MP>99?MP/100:MP)%2?1:-1)*norm);
      P2/=(((MP>99?MP/100:MP)%2?1:-1)*norm);
      P3/=norm;
      eOrGammaDyn->SetPolarization(P1,P2,P3);
      //G4cout<<eOrGammaDyn->GetPolarization()<<" ---> "<<P1*P1+P2*P2+P3*P3<<G4endl;
      //--------------------------------------------
      // The old, broken algorithm
      // G4complex tauplusplus(0,0),tauplusminus(0,0),
      // 	tauminusplus(0,0),tauminusminus(0,0);
      // const G4LevelManager *lman = G4NuclearLevelData::GetInstance()->GetLevelManager(parentZ,parentA);
      // size_t pIndex = lman->NearestLevelIndex(predecayEnergy);
      // size_t dIndex = lman->NearestLevelIndex(parentNucleus.GetExcitationEnergy()-eOrGamma->GetMomentum().getT());
      // const G4NucLevel *level = lman->GetLevel(pIndex);
      // G4double mpRatio = level->MultipolarityRatio(dIndex);
      // G4int JP1=lman->SpinTwo(pIndex);
      // G4int JP2=lman->SpinTwo(dIndex);
      // G4int MP = level->TransitionType(dIndex);
      // for(size_t k=0;k<parentPol.size();k++){
      // 	G4complex tpp(0,0),tpm(0,0),tmp(0,0),tmm(0,0);//just for kappa bit
      // 	for(size_t kappa=0;kappa<parentPol[k].size();kappa++){
      // 	  if(parentPol[k][kappa]==0.) continue;
      // 	  tpp+=parentPol[k][kappa]*G4Clebsch::WignerLittleD(2*k,2*kappa,0,gammaDir.cosTheta());
      // 	  tpm+=parentPol[k][kappa]*G4Clebsch::WignerLittleD(2*k,2*kappa,4,gammaDir.cosTheta())*G4complex(cos(2*gammaDir.phi()),sin(2*gammaDir.phi()));
      // 	  tmp+=parentPol[k][kappa]*G4Clebsch::WignerLittleD(2*k,2*kappa,-4,gammaDir.cosTheta())*G4complex(cos(2*gammaDir.phi()),-sin(2*gammaDir.phi()));
      // 	  tmm+=parentPol[k][kappa]*G4Clebsch::WignerLittleD(2*k,2*kappa,0,gammaDir.cosTheta());
      // 	}
      // 	G4PolarizationTransition polT;
      // 	double mixfactor,samefactor;
      // 	if(MP<99){//Not mixed transition
      // 	  int Lbar=MP/2;
      // 	  mixfactor = (k%2?-1:1)*sqrt(2*k+1)*(MP%2?-1:1)*G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,2,2*k)/G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,-2,2*k)*polT.FCoefficient(k,Lbar,Lbar,JP2,JP1);
      // 	  samefactor = sqrt(2*k+1)*(polT.FCoefficient(k,Lbar,Lbar,JP2,JP1));
      // 	}
      // 	else{
      // 	  int Lbar=MP/200;
      // 	  mixfactor = (k%2?-1:1)*sqrt(2*k+1)*(MP%2?-1:1)*
      // 	    (
      // 	     G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,2,2*k)/G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*Lbar,-2,2*k)*polT.FCoefficient(k,Lbar,Lbar,JP2,JP1)
      // 	     -mpRatio*(1-(k%2?-1:1))*G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*(Lbar+1),2,2*k)/G4Clebsch::ClebschGordanCoeff(2*Lbar,2,2*(Lbar+1),-2,2*k)*polT.FCoefficient(k,Lbar,Lbar+1,JP2,JP1)
      // 	     +mpRatio*mpRatio*G4Clebsch::ClebschGordanCoeff(2*(Lbar+1),2,2*(Lbar+1),2,2*k)/G4Clebsch::ClebschGordanCoeff(2*(Lbar+1),2,2*(Lbar+1),-2,2*k)*polT.FCoefficient(k,Lbar+1,Lbar+1,JP2,JP1)
      // 	     );
      // 	  samefactor = sqrt(2*k+1)*
      // 	    (
      // 	     polT.FCoefficient(k,Lbar,Lbar,JP2,JP1)
      // 	     +2*mpRatio*polT.FCoefficient(k,Lbar,Lbar+1,JP2,JP1)
      // 	     +mpRatio*mpRatio*polT.FCoefficient(k,Lbar+1,Lbar+1,JP2,JP1)
      // 	     );
      // 	}
      // 	tauplusplus+=(k%2?-1:1)*samefactor*tpp;
      // 	tauplusminus+=mixfactor*tpm;
      // 	tauminusplus+=mixfactor*tmp;
      // 	tauminusminus+=samefactor*tmm;
      // }
      // printf("\n[+][+](%lf,%lf)\n",tauplusplus.real(),tauplusplus.imag());
      // printf("[+][-](%lf,%lf)\n",tauplusminus.real(),tauplusminus.imag());
      // printf("[-][+](%lf,%lf)\n",tauminusplus.real(),tauminusplus.imag());
      // printf("[-][-](%lf,%lf)\n",tauminusminus.real(),tauminusminus.imag());
      // G4ThreeVector Stokes;
      // if((tauplusplus.imag()+tauminusminus.imag())!=0.){
      // 	G4cout<<"Complex normalization of Stokes parameters!"<<G4endl;
      // }
      // G4double norm = tauplusplus.real()+tauminusminus.real();
      // printf("norm:(%lf,%lf)\n",tauplusplus.real()+tauminusminus.real(),tauplusplus.imag()+tauminusminus.imag());
      // if((tauplusminus.imag()+tauminusplus.imag())!=0.){
      // 	G4cout<<"Stokes parameter 1 is complex!"<<G4endl;
      // }
      // Stokes.setX(-(tauplusminus.real()+tauminusplus.real())/norm);
      // printf("P1:(%lf,%lf)\n",-tauplusminus.real()-tauminusplus.real(),-tauplusminus.imag()-tauminusplus.imag());
      // if((tauminusplus.real()-tauplusminus.real())!=0.){
      // 	G4cout<<"Stokes parameter 2 is complex!"<<G4endl;
      // }
      // Stokes.setY((tauplusminus.imag()-tauminusplus.imag())/norm);
      // printf("P2:(%lf,%lf)\n",tauplusminus.imag()-tauminusplus.imag(),tauminusplus.real()-tauplusminus.real());
      // if((tauplusplus.imag()-tauminusminus.imag())!=0.){
      // 	G4cout<<"Stokes parameter 3 is complex!"<<G4endl;
      // }
      // Stokes.setZ((tauplusplus.real()-tauminusminus.real())/norm);
      // printf("P3:(%lf,%lf)\n",tauplusplus.real()-tauminusminus.real(),tauplusplus.imag()-tauminusminus.imag());
      // eOrGammaDyn->SetPolarization(Stokes.x(),Stokes.y(),Stokes.z());
      //G4cout<<Stokes<<G4endl;
      //--------------------------------------------
    }
    products->PushProducts(eOrGammaDyn);
    delete eOrGamma;

    // Now do atomic relaxation if e- is emitted
    if (applyARM) {
      G4int shellIndex = photonEvaporation->GetVacantShellNumber();
      if (shellIndex > -1) {
        G4VAtomDeexcitation* atomDeex =
          G4LossTableManager::Instance()->AtomDeexcitation();
        if (atomDeex->IsFluoActive() && parentZ > 5 && parentZ < 100) {
          G4int nShells = G4AtomicShells::GetNumberOfShells(parentZ);
          if (shellIndex >= nShells) shellIndex = nShells;
          G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIndex);
          const G4AtomicShell* shell = atomDeex->GetAtomicShell(parentZ, as);
          std::vector<G4DynamicParticle*> armProducts;

          // VI, SI
          // Allows fixing of Bugzilla 1727
          G4double deexLimit = 0.1*keV;
          if (G4EmParameters::Instance()->DeexcitationIgnoreCut())  deexLimit =0.;
          //

          atomDeex->GenerateParticles(&armProducts, shell, parentZ, deexLimit,
                                                                    deexLimit);
          G4double productEnergy = 0.;
          for (G4int i = 0; i < G4int(armProducts.size()); i++)
            productEnergy += armProducts[i]->GetKineticEnergy();

          G4double deficit = shell->BindingEnergy() - productEnergy;
          if (deficit > 0.0) { 
            // Add a dummy electron to make up extra energy
            G4double cosTh = 1.-2.*G4UniformRand();
            G4double sinTh = std::sqrt(1.- cosTh*cosTh);
            G4double phi = twopi*G4UniformRand();
         
            G4ThreeVector electronDirection(sinTh*std::sin(phi),
                                            sinTh*std::cos(phi), cosTh);
            G4DynamicParticle* extra =
              new G4DynamicParticle(G4Electron::Electron(), electronDirection,
                                    deficit);
            armProducts.push_back(extra);
          } 

          G4int nArm = armProducts.size();
          if (nArm > 0) {
            G4ThreeVector bst = dynDaughter->Get4Momentum().boostVector();
            for (G4int i = 0; i < nArm; ++i) {
              G4DynamicParticle* dp = armProducts[i];
              G4LorentzVector lv = dp->Get4Momentum().boost(bst);
              dp->Set4Momentum(lv);
              products->PushProducts(dp);
            }
          }
        }
      }
    } // if ARM on 
  } // eOrGamma

  products->PushProducts(dynDaughter);

  // Energy conservation check
  /*
  G4int newSize = products->entries();
  G4DynamicParticle* temp = 0;
  G4double KEsum = 0.0;
  for (G4int i = 0; i < newSize; i++) {
    temp = products->operator[](i);
    KEsum += temp->GetKineticEnergy();
  }
  G4double eCons = G4MT_parent->GetPDGMass() - dynDaughter->GetMass() - KEsum;
  G4cout << " IT check: Ediff (keV) = " << eCons/keV << G4endl; 
  */
  return products;
}


void G4ITDecay::DumpNuclearInfo()
{
  G4cout << " G4ITDecay for parent nucleus " << GetParentName() << G4endl;
  G4cout << " decays to " << GetDaughterName(0)
         << " + gammas (or electrons), with branching ratio " << GetBR()
         << "% and Q value " << transitionQ << G4endl;
}

#endif
