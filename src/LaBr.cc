#include "LaBr.hh"

LaBr::LaBr(Materials* mat)
{
  materials=mat;

  LaBr_Z = 3*cm; // Dirk's estimate for 2016 measurements
  
  // Dimensions based on Mirion LABR-1.5x1.5 LaBr3(Ce) Scintillation Detector
  // Spec Sheet
  
  LaBr_R  = 0.75*2.54*cm;
  LaBr_Dz = 0.75*2.54*cm;

  LaBrMaterial = materials->FindMaterial("LaBr3");

  LaBrPMT_R  = 1.156*2.54*cm;
  LaBrPMT_Dz = 2.25*2.54*cm;

  LaBrPMTMaterial = materials->FindMaterial("Al");

  // G4Polycone is defined relatve to 1 end instead of the center.
  LaBrCap_Thick = 1.0*mm;                        // wall thickness
  LaBrCap_R     = LaBr_R  + 2.*LaBrCap_Thick;    // outer
  LaBrCap_L     = 2.*LaBr_Dz + 2.*LaBrCap_Thick; // outer
  
  LaBrCapMaterial = materials->FindMaterial("Al");

  LaBrPos    = new G4ThreeVector(0., 0.,
				 -LaBr_Dz - 2.*LaBrCap_Thick - LaBr_Z);
  LaBrCapPos = new G4ThreeVector(0., 0., -LaBrCap_L - LaBr_Z);
  LaBrPMTPos = new G4ThreeVector(0., 0., 
				 -2.*LaBr_Dz-2.*LaBrCap_Thick - LaBrPMT_Dz
				 - LaBr_Z);
}

LaBr::~LaBr()
{;}
//-----------------------------------------------------------------------------
G4VPhysicalVolume* LaBr::Construct(G4LogicalVolume* experimentalHall_log)
{
  expHall_log=experimentalHall_log;
  // Crystal
  LaBr_Crys = new G4Tubs("LaBr", 0., LaBr_R, LaBr_Dz, 0., 360.*deg);
  
  LaBr_log = new G4LogicalVolume(LaBr_Crys, LaBrMaterial,
				 "LaBr_log", 0, 0, 0);

  LaBr_phys = new G4PVPlacement(G4Transform3D(NoRot,*LaBrPos),
				LaBr_log, "LaBr",
				expHall_log, false, 0);

  G4Color red (1, 0, 0); 
  G4VisAttributes* LaBrVis = new G4VisAttributes(red);
  LaBrVis->SetVisibility(true);
  LaBrVis->SetForceSolid(true);

  LaBr_log->SetVisAttributes(LaBrVis);

  // Cap
  G4double zSlice[4];
  zSlice[0] = 0;
  zSlice[1] = LaBrCap_L - LaBrCap_Thick;
  zSlice[2] = LaBrCap_L - LaBrCap_Thick;
  zSlice[3] = LaBrCap_L;

  G4double iRad[4];
  iRad[0]   = LaBrCap_R - LaBrCap_Thick;
  iRad[1]   = LaBrCap_R - LaBrCap_Thick;
  iRad[2]   = 0;
  iRad[3]   = 0;

  G4double oRad[4];
  oRad[0]   = LaBrCap_R;
  oRad[1]   = LaBrCap_R;
  oRad[2]   = LaBrCap_R;
  oRad[3]   = LaBrCap_R;

  LaBrCap = new G4Polycone("LaBrCap", 0.*deg, 360.*deg,
			   4, zSlice, iRad, oRad );
  
  LaBrCap_log = new G4LogicalVolume(LaBrCap, LaBrCapMaterial,
				    "LaBrCap_log", 0, 0, 0);

  LaBrCap_phys = new G4PVPlacement(G4Transform3D(NoRot, *LaBrCapPos),
				   LaBrCap_log, "LaBrCap",
				   expHall_log, false, 0);
  
  G4VisAttributes* CapVis = new G4VisAttributes( G4Color(1, 1, 1) );
  CapVis->SetVisibility(true);
  CapVis->SetForceWireframe(true);

  LaBrCap_log->SetVisAttributes(CapVis);

  // PMT
  LaBrPMT = new G4Tubs("LaBrPMT", 0., LaBrPMT_R, LaBrPMT_Dz,
			 0., 360.*deg);
  
  LaBrPMT_log = new G4LogicalVolume(LaBrPMT, LaBrPMTMaterial,
				 "LaBrPMT_log", 0, 0, 0);

  LaBrPMT_phys = new G4PVPlacement(G4Transform3D(NoRot,*LaBrPMTPos),
				   LaBrPMT_log, "LaBrPMT",
				   expHall_log, false, 0);

  G4Color lblue (0.0, 1.0, 1.0, 0.3); 
  G4VisAttributes* PMTVis = new G4VisAttributes(lblue);
  PMTVis->SetVisibility(true);
  PMTVis->SetForceSolid(false);

  LaBrPMT_log->SetVisAttributes(PMTVis);

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
