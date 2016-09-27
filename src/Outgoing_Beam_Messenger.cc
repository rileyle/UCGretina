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

  //REMOVE
  // LvlCmd = new G4UIcmdWithAString("/BeamOut/LevelSchemeFile",this);
  // LvlCmd->SetGuidance("Set the level scheme filename.");
  // LvlCmd->SetParameterName("choice",false);
  // LvlCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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

  TExFCmd = new G4UIcmdWithADouble("/BeamOut/TargetExcitationFraction",this);
  TExFCmd->SetGuidance("Set fraction of target excitation as compared to the projectile excitation.");
  TExFCmd->SetParameterName("choice",false);
  TExFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //REMOVE
  // tauCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/tau",this);
  // tauCmd->SetGuidance("Set lifetime of the excited state of the reaction product.");
  // tauCmd->SetParameterName("choice",false);
  // tauCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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
  
  // TB added sigma_a and sigma_b instead of just sigma
  DistSigACmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/AngDistSigmaA",this); //LR
  DistSigACmd->SetGuidance("Set sigma coefficient for Gaussian distribution of ions scattered on the target"); //LR
  DistSigACmd->SetParameterName("choice",false); //LR
  DistSigACmd->AvailableForStates(G4State_PreInit,G4State_Idle); //LR

  DistSigBCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/AngDistSigmaB",this); //LR
  DistSigBCmd->SetGuidance("Set sigma coefficient for Gaussian distribution of ions scattered on the target"); //LR
  DistSigBCmd->SetParameterName("choice",false); //LR
  DistSigBCmd->AvailableForStates(G4State_PreInit,G4State_Idle); //LR

  XsectCmd = new G4UIcmdWithAString("/BeamOut/XsectFile",this);
  XsectCmd->SetGuidance("Set the differential cross section filename.");
  XsectCmd->SetParameterName("choice",false);
  XsectCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  //REMOVE
  // QDir = new G4UIdirectory("/BeamOut/Q/");
  // QDir->SetGuidance("Charge state control for the outgoing beam.");

  // NQCmd = new G4UIcmdWithAnInteger("/BeamOut/Q/SetNumberOfChargeStates",this);
  // NQCmd->SetGuidance("Set number of charge states for the outgoing beam.");
  // NQCmd->SetParameterName("choice",false);
  // NQCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  // SQCmd = new G4UIcmdWithAnInteger("/BeamOut/Q/ChargeStateSelect",this);
  // SQCmd->SetGuidance("Select a charge state to be setup.");
  // SQCmd->SetParameterName("choice",false);
  // SQCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // SCCmd = new G4UIcmdWithAnInteger("/BeamOut/Q/Charge",this);
  // SCCmd->SetGuidance("Set charge for the charge state to be setup.");
  // SCCmd->SetParameterName("choice",false);
  // SCCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // QUFCmd = new G4UIcmdWithADouble("/BeamOut/Q/UnReactedFraction",this);
  // QUFCmd->SetGuidance("Set a fraction of the selected charge state in the unreacted beam");
  // QUFCmd->SetParameterName("choice",false);
  // QUFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  // QRFCmd = new G4UIcmdWithADouble("/BeamOut/Q/ReactedFraction",this);
  // QRFCmd->SetGuidance("Set a fraction of the selected charge state in the reacted beam");
  // QRFCmd->SetParameterName("choice",false);
  // QRFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  // QKECmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/Q/KE",this);
  // QKECmd->SetGuidance("Set kinetic energy of the central S800 trajectory for the selected charge state");
  // QKECmd->SetParameterName("choice",false);
  // QKECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  // QKEuCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/Q/KEu",this);
  // QKEuCmd->SetGuidance("Set kinetic energy per nucleon of the central S800 trajectory for the selected charge state");
  // QKEuCmd->SetParameterName("choice",false);
  // QKEuCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  UpdateCmd = new G4UIcmdWithoutParameter("/BeamOut/Update",this);
  UpdateCmd->SetGuidance("Set decay properties and report parameters for the incoming and outgoing beams.");
  
}



Outgoing_Beam_Messenger::~Outgoing_Beam_Messenger()
{
  delete BeamOutDir;
  //  delete tauCmd;
  delete ExCmd;
  delete SrcCmd;
  delete LvlCmd;
  delete TACmd;
  delete TZCmd;
  delete TExCmd;
  delete TExFCmd;
  delete RepCmd;
  delete DZCmd;
  delete DACmd;
  delete ThMinCmd;
  delete ThMTCmd;
  delete DistSigACmd;
  delete DistSigBCmd;
  delete XsectCmd;
  //REMOVE
  // delete QDir;
  // delete NQCmd;
  // delete SQCmd;
  // delete SCCmd;
  // delete QUFCmd;
  // delete QRFCmd;
  // delete QKECmd;
  // delete QKEuCmd;
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
  //REMOVE
  // if( command == LvlCmd )
  //   { BeamOut->setLvlSchemeFile(newValue);}
  if( command == LvlCmd )
    { BeamOut->setLvlDataFile(newValue);}
  if( command == TACmd )
    { BeamOut->setTarA(TACmd->GetNewIntValue(newValue));}
  if( command == TZCmd )
    { BeamOut->setTarZ(TACmd->GetNewIntValue(newValue));}
  if( command == TExFCmd )
    { BeamOut->setTFrac(TExFCmd->GetNewDoubleValue(newValue));}
  //  if( command == tauCmd )
  //    { BeamOut->settau(tauCmd->GetNewDoubleValue(newValue));}
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
  //REMOVE
  // if( command == NQCmd )
  //   { BeamOut->SetNQ(NQCmd->GetNewIntValue(newValue));}
  // if( command == SQCmd )
  //   { BeamOut->SelectQ(SQCmd->GetNewIntValue(newValue));}
  // if( command == SCCmd )
  //   { BeamOut->SetQCharge(SCCmd->GetNewIntValue(newValue));}
  // if( command == QUFCmd )
  //   { BeamOut->SetQUnReactedFraction(QUFCmd->GetNewDoubleValue(newValue));}
  // if( command == QRFCmd )
  //   { BeamOut->SetQReactedFraction(QRFCmd->GetNewDoubleValue(newValue));}
  // if( command == QKECmd )
  //   { BeamOut->SetQKE(QKECmd->GetNewDoubleValue(newValue));}
  // if( command == QKEuCmd )
  //   { BeamOut->SetQKEu(QKEuCmd->GetNewDoubleValue(newValue));}
  if( command == UpdateCmd )
    { BeamOut->Update();}
}
