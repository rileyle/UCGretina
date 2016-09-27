#ifndef Outgoing_Beam_Messenger_h
#define Outgoing_Beam_Messenger_h 1

#include "Outgoing_Beam.hh"
#include "globals.hh"
#include <vector>
using namespace std;
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"


class Outgoing_Beam_Messenger: public G4UImessenger
{
  public:
    Outgoing_Beam_Messenger(Outgoing_Beam*);
   ~Outgoing_Beam_Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Outgoing_Beam* BeamOut;    
    G4UIdirectory*             BeamOutDir;
     G4UIcmdWithAnInteger*     DACmd;
    G4UIcmdWithAnInteger*      DZCmd;
    G4UIcmdWithADoubleAndUnit* ExCmd;
    G4UIcmdWithoutParameter*   SrcCmd;
    G4UIcmdWithAString*        LvlCmd;
    G4UIcmdWithAnInteger*      TACmd;
    G4UIcmdWithAnInteger*      TZCmd;
    G4UIcmdWithADoubleAndUnit* TExCmd;
    G4UIcmdWithADouble*        TExFCmd;
    G4UIcmdWithoutParameter*   RepCmd;
    G4UIcmdWithADoubleAndUnit* ThMinCmd;
    G4UIcmdWithADoubleAndUnit* ThMTCmd;
    G4UIcmdWithADoubleAndUnit* DistSigACmd;
    G4UIcmdWithADoubleAndUnit* DistSigBCmd;
    G4UIcmdWithAString*        XsectCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};


#endif
