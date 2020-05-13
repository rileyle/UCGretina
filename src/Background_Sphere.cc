#include "Background_Sphere.hh"

Background_Sphere::Background_Sphere()
{
  BSrmin = 3.0*m;
  BSrmax = 3.4*m;
  Pos0 = new G4ThreeVector(0.,0.,0.);
  BackgroundSphereMaterial = G4Material::GetMaterial("G4_Galactic");
}

Background_Sphere::~Background_Sphere()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* Background_Sphere::Construct(G4LogicalVolume* experimentalHall_log)
{
  
  expHall_log=experimentalHall_log;
  
  BackgroundSphere = new G4Sphere("BackgroundSphere",
				  BSrmin,BSrmax,0.,360.*deg,
				  0., 180.*deg);
  
  BackgroundSphere_log = new G4LogicalVolume(BackgroundSphere,
					     BackgroundSphereMaterial,
					     "BackgroundSphere_log",0,0,0);

  BackgroundSphere_phys = new G4PVPlacement(G4Transform3D(NoRot,*Pos0),
					    BackgroundSphere_log,
					    "BackgroundSphere",expHall_log,
					    false,0);

  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0,0.3); 
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(false);
  Vis->SetForceSolid(false);

  BackgroundSphere_log->SetVisAttributes(Vis);

  return BackgroundSphere_phys;
}
//-----------------------------------------------------------------------------
void Background_Sphere::setRmin(G4double Rmin)
{
  if(Rmin<BSrmax)
    {
      BSrmin=Rmin;
      G4cout<<"----> Background Sphere inner radius set to "<<G4BestUnit(BSrmin,"Length")<< G4endl;;
      BackgroundSphere->SetInnerRadius(BSrmin);
    }
  else
    G4cout<<"----> inner radius "<<G4BestUnit(Rmin,"Length")<<" has to be smaller than the outer radius of "<<G4BestUnit(BSrmax,"Length")<<G4endl;
}
//-----------------------------------------------------------------------------
void Background_Sphere::setRmax(G4double Rmax)
{
  if(Rmax>BSrmin)
    {
      BSrmax=Rmax;
      G4cout<<"----> Background Sphere outer radius set to "<<G4BestUnit(BSrmax,"Length")<< G4endl;;
      BackgroundSphere->SetOuterRadius(BSrmax);
    }
  else
    G4cout<<"----> outer radius "<<G4BestUnit(Rmax,"Length")<<" has to be larger than the inner radius of "<<G4BestUnit(BSrmin,"Length")<<G4endl;
}
//-----------------------------------------------------------------------------
void Background_Sphere::Report()
{
     G4cout<<"----> Background Sphere material set to     "<<BackgroundSphereMaterial->GetName()<< G4endl;     
     G4cout<<"----> Background Sphere inner radius set to "<<G4BestUnit(BSrmin,"Length")<< G4endl;
     G4cout<<"----> Background Sphere outer radius set to "<<G4BestUnit(BSrmax,"Length")<< G4endl;

}
//---------------------------------------------------------------------
void Background_Sphere::setMaterial(G4String materialName)
{
  // search the material by its name 
  BackgroundSphereMaterial = G4Material::GetMaterial(materialName);  
  BackgroundSphere_log->SetMaterial(BackgroundSphereMaterial);
  G4cout<<"----> Background Sphere material set to     "<<BackgroundSphere_log->GetMaterial()->GetName()<< G4endl;                 
}
