#include "Outgoing_Beam_Messenger.hh"


Outgoing_Beam_Messenger::Outgoing_Beam_Messenger(Outgoing_Beam* BO)
:BeamOut(BO)
{ 
 
  BeamOutDir = new G4UIdirectory("/BeamOut/");
  BeamOutDir->SetGuidance("Outgoing beam control.");

  DACmd = new G4UIcmdWithAnInteger("/BeamOut/DA",this);
  DACmd->SetGuidance("Select the mass number change for the beam-like product.");
  DACmd->SetParameterName("choice",false);
  DACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DZCmd = new G4UIcmdWithAnInteger("/BeamOut/DZ",this);
  DZCmd->SetGuidance("Select the atomic number change for the beam-like product.");
  DZCmd->SetParameterName("choice",false);
  DZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ExCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/ProjectileExcitation",this);
  ExCmd->SetGuidance("Set the excitation energy for the beam-like reaction product (a level scheme file supersedes energy parameter).");
  ExCmd->SetParameterName("choice",false);
  ExCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TExCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/TargetExcitation",this);
  TExCmd->SetGuidance("Set the excitation energy of the target-like reaction product (a level scheme file supersedes energy parameter).");
  TExCmd->SetParameterName("choice",false);
  TExCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcCmd = new G4UIcmdWithoutParameter("/BeamOut/Source",this);
  SrcCmd->SetGuidance("Simulate a stationary source.");
  SrcCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  LvlCmd = new G4UIcmdWithAString("/BeamOut/LevelDataFile",this);
  LvlCmd->SetGuidance("Set the level data filename.");
  LvlCmd->SetParameterName("choice",false);
  LvlCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TACmd = new G4UIcmdWithAnInteger("/BeamOut/TargetA",this);
  TACmd->SetGuidance("Set target mass number.");
  TACmd->SetParameterName("choice",false);
  TACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TZCmd = new G4UIcmdWithAnInteger("/BeamOut/TargetZ",this);
  TZCmd->SetGuidance("Set target atomic number.");
  TZCmd->SetParameterName("choice",false);
  TZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //REMOVE
  // TExFCmd = new G4UIcmdWithADouble("/BeamOut/TargetExcitationFraction",this);
  // TExFCmd->SetGuidance("Set fraction of target excitation as compared to the projectile excitation.");
  // TExFCmd->SetParameterName("choice",false);
  // TExFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  RepCmd = new G4UIcmdWithoutParameter("/BeamOut/Report",this);
  RepCmd->SetGuidance("Report parameters for the outgoing beam.");

  ThMinCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/ThetaMin",this);
  ThMinCmd->SetGuidance("Set minimum scattering angle for the reaction product");
  ThMinCmd->SetParameterName("choice",false);
  ThMinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ThMTCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/ThetaMax",this);
  ThMTCmd->SetGuidance("Set maximum scattering angle for the reaction product");
  ThMTCmd->SetParameterName("choice",false);
  ThMTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  DistSigACmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/AngDistSigmaA",this);
  DistSigACmd->SetGuidance("Set sigma coefficient for Gaussian distribution of ions scattered on the target");
  DistSigACmd->SetParameterName("choice",false);
  DistSigACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DistSigBCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/AngDistSigmaB",this);
  DistSigBCmd->SetGuidance("Set sigma coefficient for Gaussian distribution of ions scattered on the target");
  DistSigBCmd->SetParameterName("choice",false);
  DistSigBCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  XsectCmd = new G4UIcmdWithAString("/BeamOut/XsectFile",this);
  XsectCmd->SetGuidance("Set the differential cross section filename.");
  XsectCmd->SetParameterName("choice",false);
  XsectCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/BeamOut/Update",this);
  UpdateCmd->SetGuidance("Set decay properties and report parameters for the incoming and outgoing beams.");
  
}



Outgoing_Beam_Messenger::~Outgoing_Beam_Messenger()
{
  delete BeamOutDir;
  delete ExCmd;
  delete SrcCmd;
  delete LvlCmd;
  delete TACmd;
  delete TZCmd;
  delete TExCmd;
  //  delete TExFCmd; //REMOVE
  delete RepCmd;
  delete DZCmd;
  delete DACmd;
  delete ThMinCmd;
  delete ThMTCmd;
  delete DistSigACmd;
  delete DistSigBCmd;
  delete XsectCmd;
  delete UpdateCmd;
}




void Outgoing_Beam_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == DACmd )
    { BeamOut->setDA(DACmd->GetNewIntValue(newValue));}
  if( command == DZCmd )
    { BeamOut->setDZ(DZCmd->GetNewIntValue(newValue));}
  if( command == ExCmd ){ 
    BeamOut->setEx(ExCmd->GetNewDoubleValue(newValue));
    BeamOut->setProjectileExcitation();
  }
  if( command == TExCmd ){ 
    BeamOut->setEx(TExCmd->GetNewDoubleValue(newValue));
    BeamOut->setTargetExcitation();
  }
  if( command == SrcCmd )
    { BeamOut->SetSource();}
  if( command == LvlCmd )
    { BeamOut->setLvlDataFile(newValue);}
  if( command == TACmd )
    { BeamOut->setTarA(TACmd->GetNewIntValue(newValue));}
  if( command == TZCmd )
    { BeamOut->setTarZ(TACmd->GetNewIntValue(newValue));}
  //REMOVE
  // if( command == TExFCmd )
  //   { BeamOut->setTFrac(TExFCmd->GetNewDoubleValue(newValue));}
  if( command == RepCmd )
    { BeamOut->Report();}
  if( command == ThMinCmd )
    { BeamOut->SetThetaMin(ThMinCmd->GetNewDoubleValue(newValue));}
  if( command == ThMTCmd )
    { BeamOut->SetThetaMax(ThMTCmd->GetNewDoubleValue(newValue));}
  if( command == DistSigACmd )
    { BeamOut->SetThetaSigmaA(DistSigACmd->GetNewDoubleValue(newValue));}
  if( command == DistSigBCmd )
    { BeamOut->SetThetaSigmaB(DistSigBCmd->GetNewDoubleValue(newValue));}
  if( command == XsectCmd )
    { BeamOut->setXsectFile(newValue);}
  if( command == UpdateCmd )
    { BeamOut->Update();}
}
