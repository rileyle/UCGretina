#ifndef Background_Sphere_Messenger_h
#define Background_Sphere_Messenger_h 1

#include "Background_Sphere.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

class Background_Sphere_Messenger: public G4UImessenger
{
  public:
    Background_Sphere_Messenger(Background_Sphere*);
   ~Background_Sphere_Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Background_Sphere* BackgroundSphere;
    
    G4UIdirectory*             BackgroundSphereDir;
    G4UIcmdWithAString*        MatCmd;  
    G4UIcmdWithADoubleAndUnit* RminCmd;
    G4UIcmdWithADoubleAndUnit* RmaxCmd;
    G4UIcmdWithoutParameter*   RepCmd;

};

#endif

