#include "S800.hh"

S800::S800(G4LogicalVolume* experimentalHall_log,Materials* mat)
{
  materials=mat;
  expHall_log=experimentalHall_log;

  S800_R  = 280.*mm;
  S800_Dz = 100.*mm;
  Pos0 = new G4ThreeVector(0., 0., 0.6*m);

  S800Material = materials->FindMaterial("Steel");
}

S800::~S800()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* S800::Construct()
{

  S800_Quad = new G4Tubs("S800", 3.*2.54*cm, S800_R, S800_Dz, 0., 360.*deg);
  
  S800_log = new G4LogicalVolume(S800_Quad, S800Material,
				 "S800_log", 0, 0, 0);

  S800_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
				S800_log, "S800",
				expHall_log, false, 0);

  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0, 0.3); 
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(false);

  S800_log->SetVisAttributes(Vis);

  return S800_phys;
}
//-----------------------------------------------------------------------------
void S800::setR(G4double R)
{
  S800_R=R;
  S800_Quad->SetOuterRadius(R);
  G4cout<<"----> S800 quadrupole radius set to "<<G4BestUnit(R,"Length")<< G4endl;
}
//-----------------------------------------------------------------------------
void S800::setLength(G4double length)
{
      S800_Dz=length/2.;
      G4cout<<"----> S800 length set to "<<G4BestUnit(2.*S800_Dz,"Length")<< G4endl;;
      S800_Quad->SetZHalfLength(S800_Dz);

}
//-----------------------------------------------------------------------------
void S800::Report()
{
     G4cout<<"----> S800 material set to     "<<S800Material->GetName()<< G4endl;     
     G4cout<<"----> S800 quadrupole radius set to "<<G4BestUnit(S800_R,"Length")<< G4endl;
     G4cout<<"----> S800 length set to       "<<G4BestUnit(2.*S800_Dz,"Length")<< G4endl;
}
//---------------------------------------------------------------------
void S800::setMaterial(G4String materialName)
{
  // search the material by its name 
  S800Material = materials->FindMaterial(materialName);  
  S800_log->SetMaterial(S800Material);
  G4cout<<"----> S800 material set to     "<<S800_log->GetMaterial()->GetName()<< G4endl;                 
}
