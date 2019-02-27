#include "Materials.hh"

 Materials:: Materials()
{
  // Elements
 
  G4Element* elementH  = new G4Element("Hydrogen",  "H",   1., 1.0079*g/mole);
  myElements.push_back(elementH);
  
  G4Element* elementD  = new G4Element("Deuterium", "D",   1., 2.0141*g/mole);
  myElements.push_back(elementD);
  
  G4Element* elementC  = new G4Element("Carbon",    "C",   6., 12.011*g/mole);
  myElements.push_back(elementC);

  G4Element* elementN  = new G4Element("Nitrogen",  "N",   7., 14.00674*g/mole);
  myElements.push_back(elementN);

  G4Element* elementO  = new G4Element("Oxygen",    "O",   8., 15.9994*g/mole);
  myElements.push_back(elementO);
  
  G4Element* elementMg = new G4Element("Magnesium", "Mg", 12., 24.3050*g/mole);
  myElements.push_back(elementMg);
  
  G4Element* elementAl = new G4Element("Aluminum", "Al", 13., 26.9815*g/mole);
  myElements.push_back(elementAl);
  
  G4Element* elementSi = new G4Element("Silicon",   "Si", 14., 28.0855*g/mole);
  myElements.push_back(elementSi);
  
  G4Element* elementTi = new G4Element("Titanium",  "Ti", 22., 47.90*g/mole);
  myElements.push_back(elementTi);
  
  G4Element* elementV  = new G4Element("Vanadium",  "V",  23., 50.9415*g/mole);
  myElements.push_back(elementV);

  G4Element* elementFe = new G4Element("Iron",      "Fe", 26., 55.845*g/mole);
  myElements.push_back(elementFe);

  G4Element* elementCo = new G4Element("Cobalt",    "Co", 27., 58.9332*g/mole);
  myElements.push_back(elementCo);

  G4Element* elementNi = new G4Element("Nickel",    "Ni", 28., 58.6934*g/mole);
  myElements.push_back(elementNi);

  G4Element* elementCu = new G4Element("Copper",    "Cu", 29., 63.55*g/mole);
  myElements.push_back(elementCu);

  G4Element* elementBr = new G4Element("Bromine",   "Br", 35., 79.904 *g/mole);
  myElements.push_back(elementBr);
  
  G4Element* elementMo = new G4Element("Molybdenum","Mo", 42., 95.94*g/mole);
  myElements.push_back(elementMo);

  G4Element* elementSn = new G4Element("Tin",       "Sn", 50.,118.710*g/mole);
  myElements.push_back(elementSn);

  G4Element* elementLa = new G4Element("Lanthanum", "La", 57.,138.90547*g/mole);
  myElements.push_back(elementLa);

  G4Element* elementTa = new G4Element("Tantalum",  "Ta", 73.,180.94788*g/mole);
  myElements.push_back(elementTa);

  G4Element* elementW  = new G4Element("Tungsten",  "W",  74.,183.84*g/mole);
  myElements.push_back(elementW);

  G4Element* elementPt = new G4Element("Platinum",  "Pt", 78.,195.08*g/mole);
  myElements.push_back(elementPt);

  G4Element* elementAu = new G4Element("Gold",      "Au", 79.,196.97*g/mole);
  myElements.push_back(elementAu);

  G4Element* elementPb = new G4Element("Lead",      "Pb", 82.,207.2*g/mole);
  myElements.push_back(elementPb);

  // Germanium isotopes
  G4Isotope* Ge70 = new G4Isotope("Ge70", 32, 70, 69.9242*g/mole);
  G4Isotope* Ge72 = new G4Isotope("Ge72", 32, 72, 71.9221*g/mole);
  G4Isotope* Ge73 = new G4Isotope("Ge73", 32, 73, 72.9235*g/mole);
  G4Isotope* Ge74 = new G4Isotope("Ge74", 32, 74, 73.9212*g/mole);
  G4Isotope* Ge76 = new G4Isotope("Ge76", 32, 76, 75.9214*g/mole);

  // Germanium defined via its isotopes
  G4Element* elGe = new G4Element("Germanium", "Ge", 5);
  elGe->AddIsotope(Ge70, 0.2123);
  elGe->AddIsotope(Ge72, 0.2766);
  elGe->AddIsotope(Ge73, 0.0773);
  elGe->AddIsotope(Ge74, 0.3594);
  elGe->AddIsotope(Ge76, 0.0744);
  myElements.push_back(elGe);
  
  // Materials

  G4Material* backWallMat = new G4Material("BackWallMaterial", 8.96*g/cm3, 1);
  backWallMat->AddElement(elementCu, 1.);
  myMaterials.push_back(backWallMat);

  G4Material* CD2 = new G4Material("CD2", 1.08*g/cm3, 2);
  CD2->AddElement(elementC, 1);
  CD2->AddElement(elementD, 2);
  myMaterials.push_back(CD2);

  G4Material* G10 = new G4Material("G10", 1.70*g/cm3, 4);
  G10->AddElement(elementSi, 1);
  G10->AddElement(elementO, 2);
  G10->AddElement(elementC, 3);
  G10->AddElement(elementH, 3);
  myMaterials.push_back(G10);
  
  //  G4Material* HpGe = new G4Material("HpGe", 32., 72.61*g/mole, 5.323*g/cm3);
  //  myMaterials.push_back(HpGe);
  G4Material* Ge = new G4Material("Germanium", 5.323 *g/cm3, 1);
  Ge->AddElement(elGe, 1);
  myMaterials.push_back(Ge);

  G4Material* Hevimet = new G4Material("Hevimet", 17.0*g/cm3, 3);
  Hevimet->AddElement(elementW,  0.90);
  Hevimet->AddElement(elementNi, 0.06);
  Hevimet->AddElement(elementCu, 0.04);
  myMaterials.push_back(Hevimet);

  G4Material* kapton = new G4Material("Kapton", 1.42*g/cm3, 4);
  kapton->AddElement(elementC, 0.691099);
  kapton->AddElement(elementN, 0.073298);
  kapton->AddElement(elementO, 0.209424);
  kapton->AddElement(elementH, 0.026178);
  myMaterials.push_back(kapton);

  G4Material* LaBr3 = new G4Material("LaBr3", 5.08 *g/cm3, 2);
  LaBr3->AddElement(elementLa, 0.25);
  LaBr3->AddElement(elementBr, 0.75);
  myMaterials.push_back(LaBr3);

  G4Material* preampMat = new G4Material("preampMat", 13., 26.982*g/mole,
					 1.35*g/cm3); //LR  (Air, Cu, and Al?)
  myMaterials.push_back(preampMat);

  G4Material* Steel = new G4Material("Steel", 7.86*g/cm3, 2);
  Steel->AddElement(elementFe, 98.5*perCent);
  Steel->AddElement(elementC,   1.5*perCent);
  myMaterials.push_back(Steel);
  
  G4Material* ssteel = new G4Material("ssteel", 7.7*g/cm3, 3);
  ssteel->AddElement(elementC, 0.04);
  ssteel->AddElement(elementFe, 0.88);
  ssteel->AddElement(elementCo, 0.08);
  myMaterials.push_back(ssteel);

  G4Material* Be = new G4Material("Be",      4., 9.012182*g/mole,  1.84*g/cm3);
  myMaterials.push_back(Be);
  
  G4Material* C  = new G4Material("C",       6., 12.011*g/mole,    2.15*g/cm3);
  myMaterials.push_back(C);
  
  G4Material* gC = new G4Material("glassyC", 6., 12.011*g/mole,    1.54*g/cm3);
  myMaterials.push_back(gC);
  
  G4Material* Al = new G4Material("Al",      13., 26.98154*g/mole, 2.70*g/cm3);
  myMaterials.push_back(Al);
  
  G4Material* Si = new G4Material("Si",      14., 28.0855*g/mole,  2.33*g/cm3);
  myMaterials.push_back(Si);
  
  G4Material* Fe = new G4Material("Fe",      26., 55.85*g/mole,    7.87*g/cm3);
  myMaterials.push_back(Fe);
  
  G4Material* Cu = new G4Material("Cu",      29., 63.55*g/mole,    8.96*g/cm3);
  myMaterials.push_back(Cu);
  
  G4Material* Nb = new G4Material("Nb",      41., 92.90638*g/mole, 8.57*g/cm3);
  myMaterials.push_back(Nb);
  
  G4Material* Au = new G4Material("Au",      79., 196.9*g/mole,   19.32*g/cm3);
  myMaterials.push_back(Au);
    
  G4Material* Ni = new G4Material("Ni",      28., 58.6934*g/mole,  8.908*g/cm3);
  myMaterials.push_back(Ni);
  
  G4Material* Sn = new G4Material("Sn",      50., 118.710*g/mole,  7.365*g/cm3);
  myMaterials.push_back(Sn);
  
  G4Material* Ir = new G4Material("Ir",      77., 192.217*g/mole,  22.65*g/cm3);
  myMaterials.push_back(Ir);
    
  G4Material* Pb = new G4Material("Pb",      82., 207.2*g/mole,    11.34*g/cm3);
  myMaterials.push_back(Pb);
  
  // define materials from the G4 NIST database
  NISTman = G4NistManager::Instance();

  G4Material* CsI          = NISTman->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  myMaterials.push_back(CsI);
  
  G4Material* MgO          = NISTman->FindOrBuildMaterial("G4_MAGNESIUM_OXIDE");
  myMaterials.push_back(MgO);
  
  G4Material* vacuum       = NISTman->FindOrBuildMaterial("G4_Galactic");
  myMaterials.push_back(vacuum);
  
  G4Material* polyethylene = NISTman->FindOrBuildMaterial("G4_POLYETHYLENE");
  myMaterials.push_back(polyethylene);
  
  G4Material* polypropylene= NISTman->FindOrBuildMaterial("G4_POLYPROPYLENE");
  myMaterials.push_back(polypropylene);
  
  G4Material* teflon       = NISTman->FindOrBuildMaterial("G4_TEFLON");
  myMaterials.push_back(teflon);
  
  G4Material* lH2          = NISTman->FindOrBuildMaterial("G4_lH2");
  myMaterials.push_back(lH2);
  
  G4Material* concrete     = NISTman->FindOrBuildMaterial("G4_CONCRETE");
  myMaterials.push_back(concrete);
  
  G4Material* air          = NISTman->FindOrBuildMaterial("G4_AIR");
  myMaterials.push_back(air);

  G4Material* Mg = NISTman->FindOrBuildMaterial("G4_Mg");
  myMaterials.push_back(Mg);

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
