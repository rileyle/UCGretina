#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "Stopwatch.hh"

class DetectorConstruction;
class Incoming_Beam;
class Outgoing_Beam;

// Creates master/worker user actions for Geant4 MT.
class ActionInitialization : public G4VUserActionInitialization {
public:
  ActionInitialization(DetectorConstruction* detector,
                       Incoming_Beam* beamIn,
                       Outgoing_Beam* beamOut,
                       bool enableStepping);
  ~ActionInitialization() override {delete stopwatch.Timer;}

  void BuildForMaster() const override;
  void Build() const override;

private:
  DetectorConstruction* fDetector;
  Incoming_Beam* fBeamIn;
  Outgoing_Beam* fBeamOut;
  bool fEnableStepping;

  mutable Stopwatch stopwatch;
};

#endif
