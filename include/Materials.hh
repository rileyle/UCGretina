#ifndef Materials_h
#define Materials_h 1

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"

class Materials 
{
public:

  Materials();
  ~Materials();
  
  G4Material* FindMaterial(G4String );
    
private:

  std::vector<G4Element*>  myElements;    // save pointers here to avoid

  std::vector<G4Material*> myMaterials;   // warnings of unused variables

  // Elements
  
  // G4Element* elementH;
  // G4Element* elementD;
  // G4Element* elementC;
  // G4Element* elementN;
  // G4Element* elementO;
  // G4Element* elementMg;
  // G4Element* elementAl;
  // G4Element* elementSi;
  // G4Element* elementTi;
  // G4Element* elementV;
  // G4Element* elementFe;
  // G4Element* elementCo;
  // G4Element* elementCu;
  // G4Element* elementNi;
  // G4Element* elementMo;
  // G4Element* elementSn;
  // G4Element* elementTa;
  // G4Element* elementW;
  // G4Element* elementPt;
  // G4Element* elementAu;
  // G4Element* elementPb;

  G4NistManager* NISTman;

  // Materials
  // G4Material* vacuum;
  // G4Material* HpGe;
  // G4Material* preampMat;
  // G4Material* G10;
  // G4Material* Hevimet;
  // G4Material* polyethylene;
  // G4Material* polypropylene;
  // G4Material* teflon;
  // G4Material* kapton;
  // G4Material* CD2;
  // G4Material* ssteel;
  // G4Material* CsI;
  // G4Material* MgO;
  // G4Material* Al;
  // G4Material* Fe;
  // G4Material* Cu;
  // G4Material* Nb;
  // G4Material* C;
  // G4Material* gC;
  // G4Material* Sn;
  // G4Material* Au;
  // G4Material* Pb;
  // G4Material* Ir;
  // G4Material* Si;
  // G4Material* Ni;
  // G4Material* Be;
  // G4Material* Bi;
  // G4Material* lH2;
  // G4Material* concrete;
  // G4Material* air;

};

#endif

