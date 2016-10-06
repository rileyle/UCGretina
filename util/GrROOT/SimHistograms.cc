////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////
////////////                       GrROOT
////////////
////////////          Purpose:
////////////                   To assist in the analysis of data from
////////////                 the gretina/S800 experimental setup.
////////////                          
////////////          Current Maintainers:
////////////                 Kathrin Wimmer  (wimmer@phys.s.u-tokyo.ac.jp)
////////////                 Eric Lunderberg (lunderberg@nscl.msu.edu)
////////////
////////////          Distribution:
////////////                   Please do not redistribute this software directly.
////////////                   If someone new wants a copy of this software,
////////////                 email one of the maintainers for the download link.
////////////                   This allows us to keep track of who has the software,
////////////                 and properly distribute updates and bug fixes.
////////////                 
////////////          Suggestions:
////////////                   We view the development of the software as a collaborative
////////////                 effort, and as such, welcome and appreciate any suggestions
////////////                 for bug fixes and improvements.
////////////
////////////          Disclaimer:
////////////                 This software is provided as-is, with no warranty.
////////////                 No current or future support is guaranteed for this software.
////////////
////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include <string.h>
#include <sys/time.h>


#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2S.h"
#include "TH1S.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCutG.h"
#include "TEnv.h"
#include "TKey.h"
#include "TDirectory.h"

#include "CommandLineInterface.hh"
#include "S800Calc.hh"
#include "GretinaCalc.hh"
#include "GretinaSim.hh"
#include "Mode3Calc.hh"
#include "Scaler.hh"
#include "SimHistograms.hh"
using namespace TMath;
using namespace std;

void SimHistograms::FillHistograms(GretinaCalc* gr, S800Calc* s800, Mode3Calc* m3c, GretinaSim* gs){

  //-------------------------------------------------------------------------
  //*************************************************************************
  //Fill the histograms here.
  //*************************************************************************
  //-------------------------------------------------------------------------
  for(UShort_t g = 0; g < gs->GetMult(); g++){
    Fill("egam_emitted",
	 8000,0,8000,gs->GetEmittedGamma(g)->GetEnergy());
    Fill("hposx_emitted",
	 100,-5,5,gs->GetEmittedGamma(g)->GetX());
    Fill("hposy_emitted",
	 100,-5,5,gs->GetEmittedGamma(g)->GetY());
    Fill("hposxy_emitted",
	 100,-5,5,gs->GetEmittedGamma(g)->GetX(),
	 100,-5,5,gs->GetEmittedGamma(g)->GetY());
    Fill("hposz_emitted",
	 1000,-50,50,gs->GetEmittedGamma(g)->GetZ());
    Fill("hphi_emitted",
	 360,-3.14159,3.14159,gs->GetEmittedGamma(g)->GetPhi());
    Fill("htheta_emitted",
	 180,0,3.14159,gs->GetEmittedGamma(g)->GetTheta());
    Fill("hphi_theta_emitted",
	 360,-3.14159,3.14159,gs->GetEmittedGamma(g)->GetPhi(),
	 180,0,3.14159,gs->GetEmittedGamma(g)->GetTheta());
    Fill("hbeta",
	 500,0,0.5,gs->GetEmittedGamma(g)->GetBeta());
  }
  Fill("hmult_emitted",
       10,0,10,gs->GetMult());

  // When we detected gamma rays ...
  if( gr->GetHit(0)->GetTS() == gs->GetTimeStamp() ){

    Fill("hmult_emitted_detected",
	 10,0,10, gs->GetMult(),
	 30,0,30, gr->GetMult());

    static TRACK* track = s800->GetTRACK();

    Fill("ata",
	 200,-100,100,track->GetATA());
    Fill("yta",
	 100,-50,50,track->GetYTA());
    Fill("bta",
	 200,-100,100,track->GetBTA());
    Fill("ata_vs_bta",
	 200,-100,100,track->GetATA(),
	 200,-100,100,track->GetBTA());
    Fill("dta",
	 200,-10,10,track->GetDTA());
    Fill("azita",
	 360,0,2*3.1415926,track->GetPhi());
    Fill("scatter",
	 500,0,0.5,track->GetTheta());
    Fill("ptot",
	 500,10,20,track->GetPtot());
    Fill("ppar",
	 500,10,20,track->GetPpar());
    Fill("ptra",
	 250,0,5,track->GetPtra());
    Fill("etot",
	 1000,2000,4000,track->GetEtot());


    double emitted_energy = gs->GetEmittedGamma(0)->GetEnergy();
//     double tot_energy = 0;
//     double tot_energy_dc = 0;
//     for(UShort_t g=0;g<gr->GetMult();g++){
//       tot_energy += gr->GetHit(g)->GetEnergy();
//       tot_energy_dc += gr->GetHit(g)->GetDCEnergy();
//     }
//     bool photopeak = fabs(emitted_energy - tot_energy) < 10;
//     Fill("missing_en",
// 	 5000,-1000,4000,emitted_energy - tot_energy);
//     if (photopeak){
//       Fill("egamdc_summed_photopeak",
// 	   10000,0,10000,tot_energy_dc);
//     } else {
//       Fill("egamdc_summed_compton",
// 	   10000,0,10000,tot_energy_dc);
//     }

    for(UShort_t g=0;g<gr->GetMult();g++){ // looping over gamma events
      HitCalc* hit = gr->GetHit(g);
      float energy = hit->GetEnergy();
      float energy_dc = hit->GetDCEnergy();
      float energy_dc_simple = hit->GetDCEnergy(fSett->TargetBeta());
      float energy_dc_simcheat = hit->GetDCEnergySimCheat();
      if(energy<1 || energy>6000)//useful events
	continue;

      bool photopeak = fabs(emitted_energy - energy) < 10;

      Fill("egam",
	   8000,0,8000,energy);
      Fill("egam_scatter",
	   500,0,0.5,track->GetTheta(),
	   8000,0,8000,energy);
      Fill("egamdc",
	   8000,0,8000,energy_dc);
      Fill("egamdc_scatter",
	   500,0,0.5,track->GetTheta(),
	   8000,0,8000,energy_dc);
      Fill("egamdcs",
	   8000,0,8000,energy_dc_simple);
      Fill("egamdcs_scatter",
	   500,0,0.5,track->GetTheta(),
	   8000,0,8000,energy_dc_simple);
      Fill("egamdc_simcheat",
	   8000,0,8000,energy_dc_simcheat);
      Fill("egam_emitted_diff",
	   5000,-1000,4000,emitted_energy-energy);

      Fill("thetadiff",
	   1000,-0.5,0.5,hit->GetThetaDiff());
      if (hit->GetUsedTrueFirstIP()){
	Fill("usedTrueFirst",
	     2,-0.5,1.5,1.0);
	if(photopeak){
	  Fill("usedTrueFirst_photopeak",
	       2,-0.5,1.5,1.0);
	}
	Fill("egamdc_tf",
	     8000,0,8000,energy_dc);
      } else {
	Fill("usedTrueFirst",
	     2,-0.5,1.5,0.0);
	if(photopeak){
	  Fill("usedTrueFirst_photopeak",
	       2,-0.5,1.5,0.0);
	}
	Fill("egamdc_ntf",
	     8000,0,8000,energy_dc);
	Fill("thetadiff_ntf",
	     1000,-0.5,0.5,hit->GetThetaDiff());
	Fill("angle_ntf",
	     500,0,0.5,hit->GetThetaFromTrue());
	Fill("dist_ntf",
	     1000,0,50,hit->GetDistanceFromTrue());
	Fill("perpdist_ntf",
	     1000,0,50,hit->GetPerpDistanceFromTrue());
	if (photopeak){
	  Fill("thetadiff_ntf_photopeak",
	       1000,-0.5,0.5,hit->GetThetaDiff());
	} else {
	  Fill("thetadiff_ntf_compton",
	       1000,-0.5,0.5,hit->GetThetaDiff());
	}
      }
      Fill("egamdc_thetadiff",
	   1000,-0.5,0.5,hit->GetThetaDiff(),
	   8000,0,8000,energy_dc);

      Fill("perpdist",
	   1000,0,50,hit->GetPerpDistanceFromTrue());
      if (photopeak){
	Fill("egam_photopeak",
	     8000,0,8000,energy);
	Fill("egamdc_photopeak",
	     8000,0,8000,energy_dc);
	Fill("egamdcs_photopeak",
	     8000,0,8000,energy_dc_simple);
	Fill("perpdist_photopeak",
	     1000,0,50,hit->GetPerpDistanceFromTrue());
      } else {
	Fill("egam_compton",
	     8000,0,8000,energy);
	Fill("egamdc_compton",
	     8000,0,8000,energy_dc);
	Fill("egamdcs_compton",
	     8000,0,8000,energy_dc_simple);
	Fill("perpdist_compton",
	     1000,0,50,hit->GetPerpDistanceFromTrue());
      }

      Fill("egamdc_grettheta",
	   4000,0,4000,energy_dc,
	   180,0,3.14159,hit->GetPosition().Theta());
      Fill("egamdcs_grettheta",
	   4000,0,4000,energy_dc_simple,
	   180,0,3.14159,hit->GetPosition().Theta());
      for(int i = 0; i<5; i++){
	double low = -5 + 2*i;
	double high = -5 + 2*(i+1);

	if (s800->GetTRACK()->GetDTA()>low &&
	    s800->GetTRACK()->GetDTA()<high){
	  Fill(Form("egamdc_grettheta_%d",i),
	       1000,0,4000,energy_dc,
	       180,0,3.14159,hit->GetPosition().Theta());
	}
      }
      Fill("egamdc_phi",
	   1000,0,4000,energy_dc,
	   180,-2*3.14159,0,s800->GetTRACK()->GetPhi());
      Fill(Form("egamdc_phi_d%d_c%d",
		fSett->Hole2Det(hit->GetHole()),
		hit->GetCrystal()),
	   1000,0,4000,energy_dc,
	   800,-2*3.14159,0,s800->GetTRACK()->GetPhi());

      Fill("egam_tgam",
	   8000,0,8000,energy,
	   150,0,150,hit->GetTime());

    }//gamma events



    for(UShort_t g=0;g<gr->GetMultAB();g++){ // looping over gamma events
      if(gr->GetHitAB(g)->GetEnergy()<1 || gr->GetHitAB(g)->GetEnergy()>6000)//useful events
	continue;

      HitCalc* hit = gr->GetHitAB(g);
      float energy = hit->GetEnergy();
      float energy_dc = hit->GetDCEnergy();
      Fill("egamAB",
	   8000,0,8000,energy);
      Fill("egamABdc",
	   8000,0,8000,energy_dc);
      Fill(Form("egamAB_n%d",hit->GetHitsAdded()),
	   8000,0,8000,energy);
      Fill(Form("egamABdc_n%d",hit->GetHitsAdded()),
	   8000,0,8000,energy_dc);

    }//gamma events


    Fill("hmult",
	 30,-0.5,29.5,gr->GetMult());

    for(int j=0; j<gr->GetMult(); j++){
      HitCalc* hit = gr->GetHit(j);

      Fill("hposxy",
	   600,-300,300,hit->GetPosition().X(),
	   600,-300,300,hit->GetPosition().Y());
      Fill("hposxz",
	   600,-300,300,hit->GetPosition().X(),
	   600,-300,300,hit->GetPosition().Z());
      Fill("hposyz",
	   600,-300,300,hit->GetPosition().Y(),
	   600,-300,300,hit->GetPosition().Z());
      Fill("hpos_phi",
	   360,-3.14159,3.14159,hit->GetPosition().Phi());
      Fill("hpos_theta",
	   180,0,3.14159,hit->GetPosition().Theta());
      Fill("hpos_phi_theta",
	   360,-3.14159,3.14159,hit->GetPosition().Phi(),
	   180,0,3.14159,hit->GetPosition().Theta());

      Fill("hgamma",
	   8000,0,8000,hit->GetEnergy());
      Fill(Form("hgamma_d%d_c%d",
		fSett->Hole2Det(hit->GetHole()),
		hit->GetCrystal()),
	   10000,0,20000,hit->GetEnergy());
      Fill(Form("hgammaDC_d%d_c%d",
		fSett->Hole2Det(hit->GetHole()),
		hit->GetCrystal()),
	   7000,0,7000,hit->GetDCEnergy());

      for (int k=j; k<gr->GetMult();k++){
	Fill("htimediff",
	     500,-2000,2000,hit->GetTS() - gr->GetHit(k)->GetTS());
	Fill("hegamtdiff",
	     500,-500,500,hit->GetTS() - gr->GetHit(k)->GetTS(),
	     2000,0,2000,gr->GetHit(k)->GetEnergy());
      }
    }

    for(int i=0; i<gr->GetMult(); i++){
      for(int j=i+1; j<gr->GetMult(); j++){
	double enhigh = max(gr->GetHit(i)->GetEnergy(),
			    gr->GetHit(j)->GetEnergy());
	double enlow = min(gr->GetHit(i)->GetEnergy(),
			   gr->GetHit(j)->GetEnergy());
	Fill("hegamegam_all",
	     4000,0,4000,enhigh,
	     4000,0,4000,enlow);
      }
    }

    if(gr->GetMult()==2){
      double higher = max(gr->GetHit(0)->GetEnergy(),
			  gr->GetHit(1)->GetEnergy());
      double lower = min(gr->GetHit(0)->GetEnergy(),
			 gr->GetHit(1)->GetEnergy());
      Fill("egam_egam",
	   2000,0,2000,higher,
	   2000,0,2000,lower);
      higher = max(gr->GetHit(0)->GetDCEnergy(),
		   gr->GetHit(1)->GetDCEnergy());
      lower = min(gr->GetHit(0)->GetDCEnergy(),
		  gr->GetHit(1)->GetDCEnergy());
      Fill("egam_egam_dc",
	   2000,0,2000,higher,
	   2000,0,2000,lower);
    }

    if(gr->GetMultAB()==2){
      double higher = max(gr->GetHitAB(0)->GetEnergy(),
			  gr->GetHitAB(1)->GetEnergy());
      double lower = min(gr->GetHitAB(0)->GetEnergy(),
			 gr->GetHitAB(1)->GetEnergy());
      Fill("egamAB_egamAB",
	   2000,0,2000,higher,
	   2000,0,2000,lower);
      higher = max(gr->GetHitAB(0)->GetDCEnergy(),
		   gr->GetHitAB(1)->GetDCEnergy());
      lower = min(gr->GetHitAB(0)->GetDCEnergy(),
		  gr->GetHitAB(1)->GetDCEnergy());
      Fill("egamAB_egamAB_dc",
	   2000,0,2000,higher,
	   2000,0,2000,lower);
    }

    Fill("hmult_ab",
	 30,-0.5,29.5,gr->GetMultAB());
    for(int j=0; j<gr->GetMultAB();j++){
      HitCalc* hit = gr->GetHitAB(j);
      Fill("hgamma_ab",
	   3500,0,3500,hit->GetEnergy());
    }
  }
}
