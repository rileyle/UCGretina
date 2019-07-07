#ifndef PrimaryGeneratorAction_Messenger_h
#define PrimaryGeneratorAction_Messenger_h 1

#include "PrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWith3Vector.hh"

class PrimaryGeneratorAction_Messenger: public G4UImessenger
{
public:
  PrimaryGeneratorAction_Messenger(PrimaryGeneratorAction*);
  ~PrimaryGeneratorAction_Messenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  PrimaryGeneratorAction*      PGA;    
  G4UIdirectory*               PGADir;
  G4UIdirectory*               SrcDir;
  G4UIdirectory*               RctDir;
  G4UIcmdWithAString*          SrcCmd;
  G4UIcmdWithADoubleAndUnit*   SrcECmd;
  G4UIcmdWithADoubleAndUnit*   SrcWLECmd;
  G4UIcmdWithADoubleAndUnit*   SrcWHECmd;
  G4UIcmdWithAnInteger*        SrcMultCmd;
  G4UIcmdWithADoubleAndUnit*   SrcXCmd;
  G4UIcmdWithADoubleAndUnit*   SrcYCmd;
  G4UIcmdWithADoubleAndUnit*   SrcZCmd;
  G4UIcmdWithADoubleAndUnit*   SrcRCmd;
  G4UIcmdWithADoubleAndUnit*   SrcDXCmd;
  G4UIcmdWithADoubleAndUnit*   SrcDYCmd;
  G4UIcmdWithADoubleAndUnit*   SrcSigXCmd;
  G4UIcmdWithADoubleAndUnit*   SrcSigYCmd;
  G4UIcmdWithADoubleAndUnit*   SrcCollAngCmd;
  G4UIcmdWith3Vector*          SrcCollDirCmd;
  G4UIcmdWithoutParameter*     PGASCmd;
  G4UIcmdWithoutParameter*     PGABCmd;
  G4UIcmdWithoutParameter*     SrcTFCmd;
  G4UIcmdWithoutParameter*     SrcTBCmd;
  G4UIcmdWithoutParameter*     SrcDFCmd;
  G4UIcmdWithoutParameter*     SrcDBCmd;
  G4UIcmdWithAString*          SrcThetaFileCmd;
  G4UIcmdWithoutParameter*     SrcRepCmd;
};


#endif
