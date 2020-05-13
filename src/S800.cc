#include "S800.hh"

S800::S800()
{
  S800_R      = 280.*mm; // Est.
  S800_Dz     = 100.*mm; // Est.

  S800GateValve_Dx     =   6.5*2.54*cm;  // Est.
  S800GateValve_Dy     =  12.0*2.54*cm;  // Est.
  S800GateValve_Dz     =  0.75*2.54*cm;  // Est.

  S800GateValve_yOffset = -5.5*2.54*cm;   // Est.
  S800GateValve_zOffset = (15.012 + 0.75)*2.54*cm + S800GateValve_Dz; // Est.
  S800_Offset           = S800GateValve_zOffset + S800GateValve_Dz
                          + 2.*2.54*cm + S800_Dz; // Est.
  
  S800_Pos          = new G4ThreeVector(0., 0., S800_Offset);
  S800GateValve_Pos = new G4ThreeVector(0.,
					S800GateValve_yOffset,
					S800GateValve_zOffset);
  
  S800Material = G4Material::GetMaterial("Steel");
}

S800::~S800()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* S800::Construct(G4LogicalVolume* experimentalHall_log)
{

  expHall_log=experimentalHall_log;

  S800_Quad = new G4Tubs("S800", 3.*2.54*cm, S800_R, S800_Dz, 0., 360.*deg);
  
  S800_log = new G4LogicalVolume(S800_Quad, S800Material,
				 "S800_log", 0, 0, 0);

  S800_phys = new G4PVPlacement(G4Transform3D(NoRot,
					      *S800_Pos),
				S800_log, "S800",
				expHall_log, false, 0);

  S800GateValve_Box = new G4Box("S800GateValve_Box",
				S800GateValve_Dx,
				S800GateValve_Dy,
				S800GateValve_Dz);

  S800GateValve_Hole = new G4Tubs("S800GateValve_Hole",
				  0., 3.*2.54*cm,
				  1.1*S800GateValve_Dz,
				  0., 360.*deg);

  G4ThreeVector *holePos = new G4ThreeVector(0.,
					     -S800GateValve_yOffset,
					     0.);

  S800GateValve = new G4SubtractionSolid("S800GateValve",
					 S800GateValve_Box,
					 S800GateValve_Hole,
					 G4Transform3D(NoRot,
						       *holePos));

  S800GateValve_log = new G4LogicalVolume(S800GateValve, S800Material,
					  "S800GateValve_log", 0, 0, 0);

  S800GateValve_phys = new G4PVPlacement(G4Transform3D(NoRot,
						       *S800GateValve_Pos),
					 S800GateValve_log, "S800GateValve",
					 expHall_log, false, 0);
  
  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0, 0.3); 
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(false);

  S800_log->SetVisAttributes(Vis);
  S800GateValve_log->SetVisAttributes(Vis);

  Report();
  
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
     G4cout<<"Including the S800 Quadrupole"<<G4endl; 
     G4cout<<" ----> S800 material set to     "<<S800Material->GetName()<< G4endl;     
     G4cout<<" ----> S800 quadrupole radius set to "<<G4BestUnit(S800_R,"Length")<< G4endl;
     G4cout<<" ----> S800 length set to       "<<G4BestUnit(2.*S800_Dz,"Length")<< G4endl;
}
//---------------------------------------------------------------------
void S800::setMaterial(G4String materialName)
{
  // search the material by its name 
  S800Material = G4Material::GetMaterial(materialName);  
  S800_log->SetMaterial(S800Material);
  G4cout<<"----> S800 material set to     "<<S800_log->GetMaterial()->GetName()<< G4endl;                 
}
