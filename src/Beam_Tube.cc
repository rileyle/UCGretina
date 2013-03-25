#ifndef LHTARGET
#include "Beam_Tube.hh"

Beam_Tube::Beam_Tube(G4LogicalVolume* experimentalHall_log,Materials* mat)
{
  materials=mat;
  expHall_log=experimentalHall_log;
  BTrmin = 3.0*2.54*cm - 0.058*2.54*cm;
  BTrmax = 3.0*2.54*cm;
  BTDz=44.*cm; //LR (approx target to gate valve flange)
  BTSPhi=0.*deg;
  BTDPhi=360*deg; 
  Pos0 = new G4ThreeVector(0.,0.,0.);
  BeamTubeMaterial = materials->FindMaterial("Al");

}

Beam_Tube::~Beam_Tube()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* Beam_Tube::Construct()
{

  BeamTube = new G4Tubs("BeamTube",BTrmin,BTrmax,BTDz,BTSPhi,BTDPhi);
  
  BeamTube_log = new G4LogicalVolume(BeamTube,BeamTubeMaterial,"BeamTube_log",0,0,0);

  BeamTube_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
			     BeamTube_log,"BeamTube",expHall_log,false,0);

  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0,0.3); 
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(false);

  BeamTube_log->SetVisAttributes(Vis);

  return BeamTube_phys;
}
//-----------------------------------------------------------------------------
void Beam_Tube::setRmin(G4double Rmin)
{
  if(Rmin<BTrmax)
    {
      BTrmin=Rmin;
      G4cout<<"----> Beam Tube inner radius set to "<<G4BestUnit(BTrmin,"Length")<< G4endl;;
      BeamTube->SetInnerRadius(BTrmin);
    }
  else
    G4cout<<"----> inner radius "<<G4BestUnit(Rmin,"Length")<<" has to be smaller than the outer radius of "<<G4BestUnit(BTrmax,"Length")<<G4endl;
}
//-----------------------------------------------------------------------------
void Beam_Tube::setRmax(G4double Rmax)
{
  if(Rmax>BTrmin)
    {
      BTrmax=Rmax;
      G4cout<<"----> Beam Tube outer radius set to "<<G4BestUnit(BTrmax,"Length")<< G4endl;;
      BeamTube->SetOuterRadius(BTrmax);
    }
  else
    G4cout<<"----> outer radius "<<G4BestUnit(Rmax,"Length")<<" has to be larger than the inner radius of "<<G4BestUnit(BTrmin,"Length")<<G4endl;
}
//-----------------------------------------------------------------------------
void Beam_Tube::setLength(G4double length)
{
      BTDz=length/2.;
      G4cout<<"----> Beam Tube length set to "<<G4BestUnit(2.*BTDz,"Length")<< G4endl;;
      BeamTube->SetZHalfLength(BTDz);

}
//-----------------------------------------------------------------------------
void Beam_Tube::Report()
{
     G4cout<<"----> Beam Tube material set to     "<<BeamTubeMaterial->GetName()<< G4endl;     
     G4cout<<"----> Beam Tube inner radius set to "<<G4BestUnit(BTrmin,"Length")<< G4endl;
     G4cout<<"----> Beam Tube outer radius set to "<<G4BestUnit(BTrmax,"Length")<< G4endl;
     G4cout<<"----> Beam Tube length set to       "<<G4BestUnit(2.*BTDz,"Length")<< G4endl;
 

}
//---------------------------------------------------------------------
void Beam_Tube::setMaterial(G4String materialName)
{
  // search the material by its name 
  BeamTubeMaterial = materials->FindMaterial(materialName);  
  BeamTube_log->SetMaterial(BeamTubeMaterial);
  G4cout<<"----> Beam Tube material set to     "<<BeamTube_log->GetMaterial()->GetName()<< G4endl;                 
}
#endif
