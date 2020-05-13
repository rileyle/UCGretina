#ifndef LHTARGET
#include "Greta_Chamber.hh"

Greta_Chamber::Greta_Chamber()
{

  Rmin = 178.5*mm;
  Rmax = 180.0*mm;

  BTrmin = 18.5*mm; 
  BTrmax = 20.0*mm;
  BTDz = (127.6/2*cm - Rmax)/2.;

  ChamberMaterial = G4Material::GetMaterial("Al");

  Cutaway = false;

}

Greta_Chamber::~Greta_Chamber()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* Greta_Chamber::Construct(G4LogicalVolume* experimentalHall_log)
{

  expHall_log=experimentalHall_log;
  
  G4double Phi0 = 0.*deg;
  G4double dPhi = 360.*deg;
  if(Cutaway){
    Phi0 = -90.*deg;
    dPhi = 180.*deg;
  }

  G4Sphere* sphere = new G4Sphere("sphere", Rmin, Rmax, 
				  Phi0, dPhi, 0., 180.*deg);

  G4Tubs* drill = new G4Tubs("drill", 0, BTrmax, 1.1*Rmax, 0., 360.*deg);

  G4SubtractionSolid* Chamber = new G4SubtractionSolid("Chamber", 
						       sphere, drill);

  G4RotationMatrix NoRot = G4RotationMatrix::IDENTITY;

  G4ThreeVector *Pos = new G4ThreeVector(0.,0.,0.);

  Chamber_log = new G4LogicalVolume(Chamber, ChamberMaterial, "Chamber_log",
				    0, 0, 0);

  Chamber_phys = new  G4PVPlacement(G4Transform3D(NoRot,*Pos),
				    Chamber_log, "Chamber_phys",
				    expHall_log, false, 0);

  G4Tubs* BeamTube = new G4Tubs("BeamTube", BTrmin, BTrmax, BTDz, 
				Phi0, dPhi);

  G4LogicalVolume* BeamTube_log = new G4LogicalVolume(BeamTube, 
						      ChamberMaterial,
						      "BeamTube_log",
						      0, 0, 0);

  Pos->setZ(sqrt(Rmin*Rmin - BTrmax*BTrmax)+BTDz);

  //  G4VPhysicalVolume* BeamTube1_phys = new G4PVPlacement(G4Transform3D(NoRot,
  new G4PVPlacement(G4Transform3D(NoRot, *Pos),
		    BeamTube_log,
		    "BeamTube",
		    expHall_log,
		    false,0);

  Pos->setZ(-sqrt(Rmin*Rmin - BTrmax*BTrmax)-BTDz);

  //  G4VPhysicalVolume* BeamTube2_phys = new G4PVPlacement(G4Transform3D(NoRot,
  new G4PVPlacement(G4Transform3D(NoRot, *Pos),
		    BeamTube_log,
		    "BeamTube",
		    expHall_log,
		    false,0);

  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0, 0.3); 
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(false);

  Chamber_log->SetVisAttributes(Vis);
  BeamTube_log->SetVisAttributes(Vis);

  return Chamber_phys;
}
//-----------------------------------------------------------------------------
void Greta_Chamber::Report()
{
     G4cout<<"----> Greta Chamber material set to     "<<ChamberMaterial->GetName()<< G4endl;     
     G4cout<<"----> Greta chamber inner radius set to "<<G4BestUnit(Rmin,"Length")<< G4endl;
     G4cout<<"----> Greta chamber outer radius set to "<<G4BestUnit(Rmax,"Length")<< G4endl;
     G4cout<<"----> Greta beam tube inner radius set to "<<G4BestUnit(BTrmin,"Length")<< G4endl;
     G4cout<<"----> Greta beam tube outer radius set to "<<G4BestUnit(BTrmax,"Length")<< G4endl;
     G4cout<<"----> Greta beam tube length set to       "<<G4BestUnit(2.*BTDz,"Length")<< G4endl;
     if(Cutaway)
       G4cout<<"----> Constructing the cutaway view. For visualization only!"<< G4endl;
}
#endif
