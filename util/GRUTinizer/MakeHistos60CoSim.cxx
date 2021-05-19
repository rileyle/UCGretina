
#include "TRuntimeObjects.h"

#include <iostream>
#include <map>
#include <cstdio>

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TRandom.h>
#include <TObject.h>
#include <TLine.h>

#include "TGretina.h"
#include "TS800.h"
#include "TS800.h"
#include "TGretSim.h"
#include "GValue.h"


#include "TChannel.h"

#define Q1 15
#define Q2 7
#define Q3 8
#define Q4 16
#define Q5 9
#define Q6 14
#define Q7 17
#define Q8 6
#define Q9 19

//#define BETA .37

std::map<int,int> HoleQMap;
std::map<int,std::string> LayerMap;

void InitMap() {
  HoleQMap[Q1] = 1;
  HoleQMap[Q2] = 2;
  HoleQMap[Q3] = 3;
  HoleQMap[Q4] = 4;
  HoleQMap[Q5] = 5;
  HoleQMap[Q6] = 6;
  HoleQMap[Q7] = 7;
  HoleQMap[Q8] = 8;
  HoleQMap[Q9] = 9;

  LayerMap[0] = "alpha";
  LayerMap[1] = "beta";
  LayerMap[2] = "gamma";
  LayerMap[3] = "delta";
  LayerMap[4] = "epsilon";
  LayerMap[5] = "phi";

}

#define INTEGRATION 128.0

Double_t measuredE(Double_t energy){
  Double_t resPar1 = GValue::Value("RESOLUTION_PAR_1");
  if(std::isnan(resPar1))
    resPar1 = 1.2;
  Double_t resPar2 = GValue::Value("RESOLUTION_PAR_2");
  if(std::isnan(resPar2))
    resPar2 = 0.0005;
  
  return energy + gRandom->Gaus(0.0, resPar1*sqrt(1 + resPar2*energy));
}

// extern "C" is needed to prevent name mangling.
// The function signature must be exactly as shown here,
//   or else bad things will happen.
extern "C"
void MakeHistograms(TRuntimeObjects& obj) {
  InitMap();
  TGretina *gretina = obj.GetDetector<TGretina>();
  //  TS800 *s800       = obj.GetDetector<TS800>();
  TGretSim *gretSim = obj.GetDetector<TGretSim>();
  
  Int_t    energyNChannels = 4000;
  Double_t energyLlim = 0.;
  Double_t energyUlim = 4000.;
  Double_t cosTheta  = -1;
  Double_t mCosTheta = -1;
  Double_t x1 = -1;
  Double_t y1 = -1;
  Double_t z1 = -1;
  Double_t x2 = -1;
  Double_t y2 = -1;
  Double_t z2 = -1;
  
  if(gretSim){
    for(unsigned int x=0; x<gretSim->Size(); x++){
      TGretSimHit hit = gretSim->GetGretinaSimHit(x);
      obj.FillHistogram("sim","emitted_energy",
			energyNChannels, energyLlim, energyUlim,
			hit.GetEn());
      obj.FillHistogram("sim","emitted_theta",
			180, 0., 180.,
			hit.GetTheta()*TMath::RadToDeg());
      obj.FillHistogram("sim","emitted_phi",
			360, 0., 360.,
			hit.GetPhi()*TMath::RadToDeg());
      obj.FillHistogram("sim","emitted_z",
			1000,-50., 50.,
			hit.GetZ());
      if(hit.GetEn() > 1331 && hit.GetEn() < 1334){
	x1 = sin(hit.GetTheta())*cos(hit.GetPhi());
	y1 = sin(hit.GetTheta())*sin(hit.GetPhi());
	z1 = cos(hit.GetTheta());
      }
      if(hit.GetEn() > 1171 && hit.GetEn() < 1175){
	x2 = sin(hit.GetTheta())*cos(hit.GetPhi());
	y2 = sin(hit.GetTheta())*sin(hit.GetPhi());
	z2 = cos(hit.GetTheta());
      }
    }
    cosTheta
      = (x1*x2+y1*y2+z1*z2)/(x1*x1+y1*y1+z1*z1)/(x2*x2+y2*y2+z2*z2);
    obj.FillHistogram("sim","emitted_delta",
		      100, -1, 1.,
		      cosTheta);
  }

  if(!gretina)
    return;
  
  Double_t calorimeterEnergy = 0.;
  std::vector<TGretinaHit> hits;

  // Gamma-gated crystal spectrum
  Double_t eGateLlim = 1327.5;
  Double_t eGateUlim = 1337.5;
  int iGate = -1;
  for(unsigned int i=0; i<gretina->Size(); i++){
    TGretinaHit hit = gretina->GetGretinaHit(i);
    if( hit.GetCoreEnergy() > eGateLlim &&
        hit.GetCoreEnergy() < eGateUlim )
      iGate = i;
  }
  if(iGate>=0){
    for(unsigned int i=0; i<gretina->Size(); i++){
      TGretinaHit hit = gretina->GetGretinaHit(i);
      if(i != iGate){
	Double_t mE = measuredE(hit.GetCoreEnergy());
	obj.FillHistogram("energy",
			  Form("energy_%.0f",
			       (eGateLlim+eGateUlim)/2.0),
			  energyNChannels, energyLlim, energyUlim, mE);
	obj.FillHistogram("energy",
			  Form("fold_vs_energy_%.0f",
			       (eGateLlim+eGateUlim)/2.0),
			  energyNChannels/8, energyLlim, energyUlim, mE,
			  20, 0, 20, hit.NumberOfInteractions());

	// Count segment fold
	int segment_fold = hit.NumberOfInteractions();
	for(int y=0; y < hit.NumberOfInteractions(); y++)
	  for(int z = y+1; z < hit.NumberOfInteractions(); z++)
	    if(hit.GetSegmentId(y) == hit.GetSegmentId(z)){
	      segment_fold--;
	      break;
	    }

	obj.FillHistogram("energy",
			  Form("segfold_vs_energy_%.0f",
			       (eGateLlim+eGateUlim)/2.0),
			  energyNChannels/8, energyLlim, energyUlim, mE,
			  20, 0, 20, segment_fold);
	
      }
    }
  }

  // Addback preprocessing
  for(unsigned int x=0; x<gretina->Size(); x++){

    TGretinaHit hit = gretina->GetGretinaHit(x);
    
    if(hit.GetCoreEnergy() > energyLlim &&
       hit.GetCoreEnergy() < energyUlim){

      calorimeterEnergy += measuredE(hit.GetCoreEnergy());
      hits.push_back(hit);

    }
  }
  
  int max_layer = -1;

  for(unsigned int x=0; x<gretina->Size(); x++){

    TGretinaHit hit = gretina->GetGretinaHit(x);
    Double_t mE = measuredE(hit.GetCoreEnergy());

    //                directory, histogram
    obj.FillHistogram("energy", "overview",
		      energyNChannels, energyLlim, energyUlim, mE,
		      100, 0, 100, hit.GetCrystalId());

    obj.FillHistogram("energy", "energy",
		      energyNChannels, energyLlim, energyUlim, mE);

    obj.FillHistogram("energy", "fold_vs_energy",
		      energyNChannels/8, energyLlim, energyUlim, mE,
		      20, 0, 20, hit.NumberOfInteractions());

    // Symmetrized gamma-gamma matrix and angular correlation
    for(unsigned int y=x+1; y<gretina->Size(); y++){
      TGretinaHit hit2 = gretina->GetGretinaHit(y);
      Double_t mE2 = measuredE(hit2.GetCoreEnergy());
      obj.FillHistogram("energy", "gamma_gamma",
			energyNChannels/4, energyLlim, energyUlim, mE,
			energyNChannels/4, energyLlim, energyUlim, mE2);
      obj.FillHistogram("energy", "gamma_gamma",
			energyNChannels/4, energyLlim, energyUlim, mE2,
			energyNChannels/4, energyLlim, energyUlim, mE);

      if( (mE > 1168. && mE < 1178. && mE2 > 1328. && mE2 < 1337.)
	  || (mE2 > 1168. && mE2 < 1178. && mE > 1328. && mE < 1337.) ){
	//	TVector3 r1 = hit.GetFirstIntPosition_2();
	//	TVector3 r2 = hit2.GetFirstIntPosition_2();
	TVector3 r1 = hit.GetPosition();
	TVector3 r2 = hit2.GetPosition();
        mCosTheta = r1.Dot(r2)/r1.Mag()/r2.Mag();
	obj.FillHistogram("position","delta",
			  100, -1, 1.,
			  mCosTheta);
	obj.FillHistogram("position","cosTheta_mCosTheta",
			  100, -1, 1.,
			  cosTheta,
			  100, -1, 1.,
			  mCosTheta);
      }

    }
    
    // Count segment fold
    int segment_fold = hit.NumberOfInteractions();
    for(int y=0; y < hit.NumberOfInteractions(); y++)
      for(int z = y+1; z < hit.NumberOfInteractions(); z++)
	if(hit.GetSegmentId(y) == hit.GetSegmentId(z)){
	  segment_fold--;
	  break;
	}

    obj.FillHistogram("energy",
		      "segfold_vs_energy",
		      energyNChannels/8, energyLlim, energyUlim, mE,
		      20, 0, 20, segment_fold);

    if( hit.GetCrystalId()%2 )
      obj.FillHistogram("energy",  "energy_A",
			energyNChannels, energyLlim, energyUlim, mE);
    else
      obj.FillHistogram("energy",  "energy_B", 
			energyNChannels, energyLlim, energyUlim, mE);
                        
    // Peter wrote these to give Dirk his Phi range (0-360).
    obj.FillHistogram("position", "theta_vs_phi",
		      360, 0., 360.,
		      hit.GetPhi()*TMath::RadToDeg(),
		      180, 0., 180.,
		      hit.GetTheta()*TMath::RadToDeg());

    // Position spectra in crystal coordinates
    if(hit.NumberOfInteractions()){
      obj.FillHistogram("position",
			Form("crys_%d_x", hit.GetCrystalId()),
			1200, -60, 60, hit.GetLocalIntPosition(0).X());

      obj.FillHistogram("position",
			Form("crys_%d_y", hit.GetCrystalId()),
			1200, -60, 60, hit.GetLocalIntPosition(0).Y());

      obj.FillHistogram("position",
			Form("crys_%d_z", hit.GetCrystalId()),
			1200, -10, 100, hit.GetLocalIntPosition(0).Z());

      obj.FillHistogram("position",
			Form("crys_%d_xy", hit.GetCrystalId()),
			1200, -60, 60, hit.GetLocalIntPosition(0).X(),
			1200, -60, 60, hit.GetLocalIntPosition(0).Y());

      obj.FillHistogram("position",
			Form("crys_%d_xz", hit.GetCrystalId()),
			1200, -60, 60,  hit.GetLocalIntPosition(0).X(),
			1200, -10, 100, hit.GetLocalIntPosition(0).Z());

      obj.FillHistogram("position",
			Form("crys_%d_yz", hit.GetCrystalId()),
			1200, -60, 60,  hit.GetLocalIntPosition(0).Y(),
			1200, -10, 100, hit.GetLocalIntPosition(0).Z());
    }
    
    for(int y=0; y < hit.NumberOfInteractions(); y++){

      int layer = hit.GetSegmentId(y)/6;

      if(layer > max_layer) max_layer = layer;
      
      obj.FillHistogram("layers",
			Form("theta_vs_phi_%s", LayerMap[layer].c_str()),
			360, 0., 360.,
			hit.GetPhi()*TMath::RadToDeg(),
			180, 0., 180.,
			hit.GetTheta()*TMath::RadToDeg());
    }

    obj.FillHistogram("layers", "max_layer", 12, -2, 10,
		      max_layer);

    if(max_layer == 5){
	obj.FillHistogram("energy", "energy_involves_phi",
			  energyNChannels, energyLlim, energyUlim, mE);
	obj.FillHistogram("energy", "overview_involves_phi",
			  energyNChannels, energyLlim, energyUlim, mE,
			  100, 0, 100, hit.GetCrystalId());
    }
    
    for(int k = 5; k > 0; k--){
      if(max_layer < k){
	obj.FillHistogram("energy",
			  Form("energy_below_%s", LayerMap[k].c_str()),
			  energyNChannels, energyLlim, energyUlim, mE);
	if( hit.GetCrystalId()%2 )
	  obj.FillHistogram("energy",
			    Form("energy_A_below_%s", LayerMap[k].c_str()),
			    energyNChannels, energyLlim, energyUlim, mE);
	else
	  obj.FillHistogram("energy",
			    Form("energy_B_below_%s", LayerMap[k].c_str()),
			    energyNChannels, energyLlim, energyUlim, mE);
      }

    }
    
  }
  
  // Addback
  obj.FillHistogram("addback",  "calorimeter",
		    energyNChannels, energyLlim, energyUlim,
		    calorimeterEnergy);
  for(int k = 5; k > 0; k--){
    if(max_layer < k){
      obj.FillHistogram("addback",
			Form("calorimeter_below_%s",
			     LayerMap[k].c_str()),
			energyNChannels, energyLlim, energyUlim,
			calorimeterEnergy);
    }
  }
  
  // For energy-gated addback spectra
  std::vector<Double_t> ABenergy;
  std::vector<TString> ABtype;
  Int_t    Naddback = 0;
  iGate = -1;
  
  while(hits.size() > 0){
    TGretinaHit currentHit = hits.back();
    hits.pop_back();
    
    // Find and add all hits in a cluster of adjacent crystals including
    // the current hit.
    //
    // CAUTION: This clustering includes neighbors of neighbors!
    std::vector<TGretinaHit> cluster;
    cluster.push_back(currentHit);
    unsigned int lastClusterSize = 0;
    while(lastClusterSize < cluster.size()){
      for(unsigned int i = 0; i < cluster.size(); i++){
	for(unsigned int j = 0; j < hits.size(); j++){
	  TVector3 distance = cluster[i].GetCrystalPosition()
	                       - hits[j].GetCrystalPosition();

	  obj.FillHistogram("position",  "crystal_separation",
			    1000, 0., 1000.,
			    distance.Mag());

	  if(distance.Mag() < 80.){ // Neighbors
	    cluster.push_back(hits.back());
	    hits.pop_back();
	  }
	}
      }
      lastClusterSize = cluster.size();
    }
    
    // Calculate the total energy deposited in the cluster,
    // and count the pairs of neighbors.
    Int_t neighbors = 0;
    Double_t addbackEnergy = 0.;
    for(unsigned int i = 0; i < cluster.size(); i++){
      addbackEnergy += measuredE(cluster[i].GetCoreEnergy());
      for(unsigned int j = i+1; j < cluster.size(); j++){
	TVector3 distance =   cluster[i].GetCrystalPosition()
	                    - cluster[j].GetCrystalPosition();
	if(distance.Mag() < 80.) neighbors++;
      }
    }

    TString addbackType;
    if(neighbors == 0 && cluster.size() == 1)
      addbackType = "addback_n0";
    else if(neighbors == 1 && cluster.size() == 2)
      addbackType = "addback_n1";
    else if(neighbors == 3 && cluster.size() == 3)
      addbackType = "addback_n2";
    else
      addbackType = "addback_ng";

    // Fill addback histograms.

    obj.FillHistogram("addback",  addbackType,
		      energyNChannels, energyLlim, energyUlim,
		      addbackEnergy);

    if(addbackType == "addback_n0"
       || addbackType == "addback_n1"){
      obj.FillHistogram("addback",  "addback_n0n1",
			energyNChannels, energyLlim, energyUlim,
			addbackEnergy);
    }
    
    if(addbackType == "addback_n0"
       || addbackType == "addback_n1"
       || addbackType == "addback_n2"){
      obj.FillHistogram("addback",  "addback_n0n1n2",
			energyNChannels, energyLlim, energyUlim,
			addbackEnergy);
    }
    
    obj.FillHistogram("addback",  "clusterSize_vs_neighborPairs",
		      20, 0, 20, neighbors,
		      10, 0, 10, cluster.size());

    // For energy-gated addback spectra
    if(addbackEnergy > eGateLlim &&
       addbackEnergy < eGateUlim)
      iGate = Naddback;

    ABenergy.push_back(addbackEnergy);
    ABtype.push_back(addbackType);
    Naddback++;

  }

  // Fill energy-gated addback spectra
  if(iGate >= 0){
    for(int i = 0; i<Naddback; i++){
      if(i != iGate){
	obj.FillHistogram("addback",
			  ABtype[i]+Form("_%.0f",
					   (eGateLlim+eGateUlim)/2.0),
			  energyNChannels, energyLlim, energyUlim,
			  ABenergy[i]);

	if(ABtype[i] == "addback_n0"
	   || ABtype[i] == "addback_n1"){
	  obj.FillHistogram("addback",
			    Form("addback_n0n1_%.0f",
				 (eGateLlim+eGateUlim)/2.0),
			    energyNChannels, energyLlim, energyUlim,
			    ABenergy[i]);
	}
    
	if(ABtype[i] == "addback_n0"
	   || ABtype[i] == "addback_n1"
	   || ABtype[i] == "addback_n2"){
	  obj.FillHistogram("addback",
			    Form("addback_n0n1n2_%.0f",
				 (eGateLlim+eGateUlim)/2.0),
			    energyNChannels, energyLlim, energyUlim,
			    ABenergy[i]);
	}
	
      }
    }
  }
  
  TList *list = &(obj.GetObjects());
  int numobj = list->GetSize();
  
  if(numobj!=list->GetSize())
    list->Sort();
}
