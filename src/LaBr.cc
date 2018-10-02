#include "LaBr.hh"

LaBr::LaBr(G4LogicalVolume* experimentalHall_log,Materials* mat)
{
  materials=mat;
  expHall_log=experimentalHall_log;

  LaBr_R  = 0.75*2.54*cm;
  LaBr_Dz = 0.75*2.54*cm;

  LaBrPos = new G4ThreeVector(0., 0., -1.5*2.54*cm);

  LaBrMaterial = materials->FindMaterial("LaBr3");

  LaBrPMT_R  = 1.*2.54*cm;
  LaBrPMT_Dz = 3.*2.54*cm;

  LaBrPMTMaterial = materials->FindMaterial("Al");

  LaBrPMTPos = new G4ThreeVector(0., 0.,
				 -LaBr_Dz-LaBrPMT_Dz-1.5*2.54*cm);
  
}

LaBr::~LaBr()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* LaBr::Construct()
{


  LaBr_Crys = new G4Tubs("LaBr", 0., LaBr_R, LaBr_Dz, 0., 360.*deg);
  
  LaBr_log = new G4LogicalVolume(LaBr_Crys, LaBrMaterial,
				 "LaBr_log", 0, 0, 0);

  LaBr_phys = new G4PVPlacement(G4Transform3D(NoRot,*LaBrPos),
				LaBr_log, "LaBr",
				expHall_log, false, 0);

  LaBrPMT = new G4Tubs("LaBrPMT", 0., LaBrPMT_R, LaBrPMT_Dz,
			 0., 360.*deg);
  
  LaBrPMT_log = new G4LogicalVolume(LaBrPMT, LaBrPMTMaterial,
				 "LaBrPMT_log", 0, 0, 0);

  LaBrPMT_phys = new G4PVPlacement(G4Transform3D(NoRot,*LaBrPMTPos),
				   LaBrPMT_log, "LaBrPMT",
				   expHall_log, false, 0);

  // Visualization Attributes

  G4Colour lblue (0.0, 1.0, 1.0, 0.3); 
  G4VisAttributes* Vis = new G4VisAttributes(lblue);
  Vis->SetVisibility(true);
  Vis->SetForceSolid(false);

  LaBr_log->SetVisAttributes(Vis);
  LaBrPMT_log->SetVisAttributes(Vis);

  Report();
  
  return LaBr_phys;
}
//-----------------------------------------------------------------------------
void LaBr::MakeSensitive(TrackerGammaSD* TrackerGamma)
{
  LaBr_log->SetSensitiveDetector(TrackerGamma);
}
//---------------------------------------------------------------------
void LaBr::setR(G4double R)
{
  LaBr_R=R;
  LaBr_Crys->SetOuterRadius(R);
  G4cout<<"----> LaBr quadrupole radius set to "<<G4BestUnit(R,"Length")<< G4endl;
}
//-----------------------------------------------------------------------------
void LaBr::setLength(G4double length)
{
      LaBr_Dz=length/2.;
      G4cout<<"----> LaBr length set to "<<G4BestUnit(2.*LaBr_Dz,"Length")<< G4endl;;
      LaBr_Crys->SetZHalfLength(LaBr_Dz);

}
//-----------------------------------------------------------------------------
void LaBr::Report()
{
     G4cout<<"Including the LaBr"<<G4endl; 
     G4cout<<" ----> LaBr material set to     "<<LaBrMaterial->GetName()<< G4endl;     
     G4cout<<" ----> LaBr quadrupole radius set to "<<G4BestUnit(LaBr_R,"Length")<< G4endl;
     G4cout<<" ----> LaBr length set to       "<<G4BestUnit(2.*LaBr_Dz,"Length")<< G4endl;
}
//---------------------------------------------------------------------
void LaBr::setMaterial(G4String materialName)
{
  // search the material by its name 
  LaBrMaterial = materials->FindMaterial(materialName);  
  LaBr_log->SetMaterial(LaBrMaterial);
  G4cout<<"----> LaBr material set to     "<<LaBr_log->GetMaterial()->GetName()<< G4endl;                 
}
