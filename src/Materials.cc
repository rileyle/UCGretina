#include "Materials.hh"

 Materials:: Materials()
{
  // Elements
 
  elementH  = new G4Element("Hydrogen",  "H",  1.,  1.0079*g/mole);
  elementD  = new G4Element("Deuterium", "D",  1.,  2.0141*g/mole);
  elementC  = new G4Element("Carbon",    "C",  6.,  12.011*g/mole);
  elementN  = new G4Element("Nitrogen",  "N",  7.,  14.007*g/mole);
  elementO  = new G4Element("Oxygen",    "O",  8., 15.9994*g/mole);
  elementMg = new G4Element("Magnesium", "Mg",12., 24.3050*g/mole);
  elementAl = new G4Element("Aluminium", "Al",13., 26.9815*g/mole);
  elementSi = new G4Element("Silicon",   "Si",14., 28.0855*g/mole);
  elementTi = new G4Element("Titanium",  "Ti",22.,   47.90*g/mole);
  elementV  = new G4Element("Vanadium",  "V", 23., 50.9415*g/mole);
  elementFe = new G4Element("Iron",      "Fe",26.,  55.845*g/mole);
  elementCo = new G4Element("Cobalt",    "Co",27., 58.9332*g/mole);
  elementCu = new G4Element("Copper",    "Cu",29.,   63.55*g/mole);
  elementMo = new G4Element("Molybdenum","Mo",42.,   95.94*g/mole);
  elementSn = new G4Element("Tin",       "Sn",50., 118.710*g/mole);
  elementPt = new G4Element("Platinum",  "Pt",78.,  195.08*g/mole);
  elementAu = new G4Element("Gold",      "Au",79.,  196.97*g/mole);

  // Materials

  HpGe = new G4Material("HpGe", 32., 72.61*g/mole, 5.323*g/cm3);
  preampMat = new G4Material("preampMat", 13., 26.982*g/mole, 1.35*g/cm3);   //LR  (Air, copper, and aluminum?)

  G10 = new G4Material("G10", 1.70*g/cm3, 4);
  G10->AddElement(elementSi, 1);
  G10->AddElement(elementO, 2);
  G10->AddElement(elementC, 3);
  G10->AddElement(elementH, 3);

  CD2 = new G4Material("CD2", 1.08*g/cm3, 2);
  CD2->AddElement(elementC, 1);
  CD2->AddElement(elementD, 2);

  ssteel = new G4Material("ssteel", 7.7*g/cm3, 3);
  ssteel->AddElement(elementC, 0.04);
  ssteel->AddElement(elementFe, 0.88);
  ssteel->AddElement(elementCo, 0.08);

  kapton = new G4Material("Kapton", 1.42*g/cm3, 4);
  kapton->AddElement(elementC, 0.691099);
  kapton->AddElement(elementN, 0.073298);
  kapton->AddElement(elementO, 0.209424);
  kapton->AddElement(elementH, 0.026178);

  // define materials from the G4 NIST database
  NISTman = G4NistManager::Instance();

  CsI           = NISTman->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  MgO           = NISTman->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE");
  vacuum        = NISTman->FindOrBuildMaterial("G4_Galactic");
  polyethylene  = NISTman->FindOrBuildMaterial("G4_POLYETHYLENE");
  polypropylene = NISTman->FindOrBuildMaterial("G4_POLYPROPYLENE");
  lH2           = NISTman->FindOrBuildMaterial("G4_lH2");
  concrete      = NISTman->FindOrBuildMaterial("G4_CONCRETE");
  air           = NISTman->FindOrBuildMaterial("G4_AIR");

  Be = new G4Material("Be",       4., 9.012182*g/mole,  1.84*g/cm3);
  C  = new G4Material("C",        6., 12.011*g/mole,    2.15*g/cm3);
  gC = new G4Material("glassyC",  6., 12.011*g/mole,  1.54*g/cm3);
  Al = new G4Material("Al",       13., 26.98153*g/mole,  2.70*g/cm3);
  Si = new G4Material("Si",       14., 28.0855*g/mole,   2.33*g/cm3);
  Fe = new G4Material("Fe",       26., 55.85*g/mole,     7.87*g/cm3);
  Cu = new G4Material("Cu",       29., 63.55*g/mole,     8.96*g/cm3);
  Nb = new G4Material("Nb",       41., 92.90638*g/mole,  8.57*g/cm3);
  Au = new G4Material("Au",       79., 196.9*g/mole,    19.32*g/cm3);
  Be = new G4Material("Be",        4., 9.012182*g/mole,  1.84*g/cm3);
  Si = new G4Material("Si",       14., 28.0855*g/mole,   2.33*g/cm3);
  Sn = new G4Material("Sn",       50., 118.710*g/mole,  7.365*g/cm3);
  Ir = new G4Material("Ir",       77., 192.217*g/mole,  22.65*g/cm3);
  Au = new G4Material("Au",       79., 196.9*g/mole,    19.32*g/cm3);

}

 Materials::~ Materials()
{;}
//-----------------------------------------------------------------------------
G4Material*  Materials::FindMaterial(G4String materialName)
{

   // search the material by its name 
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);  

  return pttoMaterial;
  
}
