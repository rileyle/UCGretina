#include "Outgoing_Beam_Messenger.hh"


Outgoing_Beam_Messenger::Outgoing_Beam_Messenger(Outgoing_Beam* BO)
:BeamOut(BO)
{ 
 
  BeamOutDir = new G4UIdirectory("/BeamOut/");
  BeamOutDir->SetGuidance("Outgoing beam control.");

  DACmd = new G4UIcmdWithAnInteger("/BeamOut/DA",this);
  DACmd->SetGuidance("Select the mass number change for the outgoing beam.");
  DACmd->SetParameterName("choice",false);
  DACmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DZCmd = new G4UIcmdWithAnInteger("/BeamOut/DZ",this);
  DZCmd->SetGuidance("Select the atomic number change for the outgoing beam.");
  DZCmd->SetParameterName("choice",false);
  DZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ExCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/ProjectileExcitation",this);
  ExCmd->SetGuidance("Set excitation energy for the outgoing beam.");
  ExCmd->SetParameterName("choice",false);
  ExCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcCmd = new G4UIcmdWithoutParameter("/BeamOut/Source",this);
  SrcCmd->SetGuidance("Simulate a stationary source.");
  SrcCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  LvlCmd = new G4UIcmdWithAString("/BeamOut/LevelSchemeFile",this);
  LvlCmd->SetGuidance("Set the level scheme filename.");
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

  TExCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/TargetExcitation",this);
  TExCmd->SetGuidance("Set target excitation energy.");
  TExCmd->SetParameterName("choice",false);
  TExCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TExFCmd = new G4UIcmdWithADouble("/BeamOut/TargetExcitationFraction",this);
  TExFCmd->SetGuidance("Set fraction of target excitation as compared to the projectile excitation.");
  TExFCmd->SetParameterName("choice",false);
  TExFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  tauCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/tau",this);
  tauCmd->SetGuidance("Set lifetime of the excited state of the outgoing beam.");
  tauCmd->SetParameterName("choice",false);
  tauCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
/*
  betaCmd = new G4UIcmdWithADouble("/BeamOut/beta",this);
  betaCmd->SetGuidance("Set beta for the outgoing beam.");
  betaCmd->SetParameterName("choice",false);
  betaCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DoppZCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/DoppZ",this);
  DoppZCmd->SetGuidance("Set Z position for the outgoing beam Doppler corrections.");
  DoppZCmd->SetParameterName("choice",false);
  DoppZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DoppYCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/DoppY",this);
  DoppYCmd->SetGuidance("Set Y position for the outgoing beam Doppler corrections.");
  DoppYCmd->SetParameterName("choice",false);
  DoppYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DoppXCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/DoppX",this);
  DoppXCmd->SetGuidance("Set X position for the outgoing beam Doppler corrections.");
  DoppXCmd->SetParameterName("choice",false);
  DoppXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
*/
  RepCmd = new G4UIcmdWithoutParameter("/BeamOut/Report",this);
  RepCmd->SetGuidance("Report parameters for the outgoing beam.");

  ThMTCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/ThetaMax",this);
  ThMTCmd->SetGuidance("Set maximum scattering angle for ions Coulomb scattered on the target");
  ThMTCmd->SetParameterName("choice",false);
  ThMTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  AphTCmd = new G4UIcmdWithADouble("/BeamOut/AngDistAlpha",this);
  AphTCmd->SetGuidance("Set alpha coefficient for angular distribution of ions Coulomb scattered on the target");
  AphTCmd->SetParameterName("choice",false);
  AphTCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
    // TB added sigma_a and sigma_b instead of just sigma
  DistSigACmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/AngDistSigmaA",this); //LR
  DistSigACmd->SetGuidance("Set sigma coefficient for Gaussian distribution of ions scattered on the target"); //LR
  DistSigACmd->SetParameterName("choice",false); //LR
  DistSigACmd->AvailableForStates(G4State_PreInit,G4State_Idle); //LR

  DistSigBCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/AngDistSigmaB",this); //LR
  DistSigBCmd->SetGuidance("Set sigma coefficient for Gaussian distribution of ions scattered on the target"); //LR
  DistSigBCmd->SetParameterName("choice",false); //LR
  DistSigBCmd->AvailableForStates(G4State_PreInit,G4State_Idle); //LR
  
  QDir = new G4UIdirectory("/BeamOut/Q/");
  QDir->SetGuidance("Charge state control for the outgoing beam.");

  NQCmd = new G4UIcmdWithAnInteger("/BeamOut/Q/SetNumberOfChargeStates",this);
  NQCmd->SetGuidance("Set number of charge states for the outgoing beam.");
  NQCmd->SetParameterName("choice",false);
  NQCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SQCmd = new G4UIcmdWithAnInteger("/BeamOut/Q/ChargeStateSelect",this);
  SQCmd->SetGuidance("Select a charge state to be setup.");
  SQCmd->SetParameterName("choice",false);
  SQCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SCCmd = new G4UIcmdWithAnInteger("/BeamOut/Q/Charge",this);
  SCCmd->SetGuidance("Set charge for the charge state to be setup.");
  SCCmd->SetParameterName("choice",false);
  SCCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  QUFCmd = new G4UIcmdWithADouble("/BeamOut/Q/UnReactedFraction",this);
  QUFCmd->SetGuidance("Set a fraction of the selected charge state in the unreacted beam");
  QUFCmd->SetParameterName("choice",false);
  QUFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  QRFCmd = new G4UIcmdWithADouble("/BeamOut/Q/ReactedFraction",this);
  QRFCmd->SetGuidance("Set a fraction of the selected charge state in the reacted beam");
  QRFCmd->SetParameterName("choice",false);
  QRFCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  QKECmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/Q/KE",this);
  QKECmd->SetGuidance("Set kinetic energy of the central S800 trajectory for the selected charge state");
  QKECmd->SetParameterName("choice",false);
  QKECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  QKEuCmd = new G4UIcmdWithADoubleAndUnit("/BeamOut/Q/KEu",this);
  QKEuCmd->SetGuidance("Set kinetic energy per nucleon of the central S800 trajectory for the selected charge state");
  QKEuCmd->SetParameterName("choice",false);
  QKEuCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

	a0cmd = new G4UIcmdWithADouble("/BeamOut/seta0",this);
	a2cmd = new G4UIcmdWithADouble("/BeamOut/seta2",this);
	a4cmd = new G4UIcmdWithADouble("/BeamOut/seta4",this);
	a0Targetcmd = new G4UIcmdWithADouble("/BeamOut/setTargeta0",this);
	a2Targetcmd = new G4UIcmdWithADouble("/BeamOut/setTargeta2",this);
	a4Targetcmd = new G4UIcmdWithADouble("/BeamOut/setTargeta4",this);
}



Outgoing_Beam_Messenger::~Outgoing_Beam_Messenger()
{
  delete BeamOutDir;
//  delete DoppXCmd;
//  delete DoppYCmd;
//  delete DoppZCmd;
//  delete betaCmd;
  delete tauCmd;
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
  delete ThMTCmd;
  delete AphTCmd;
  delete DistSigACmd;
  delete DistSigBCmd;
  delete QDir;
  delete NQCmd;
  delete SQCmd;
  delete SCCmd;
  delete QUFCmd;
  delete QRFCmd;
  delete QKECmd;
  delete QKEuCmd;
	delete a0Targetcmd;
	delete a2Targetcmd;
	delete a4Targetcmd;
	delete a0cmd;
	delete a2cmd;
	delete a4cmd;
}




void Outgoing_Beam_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 

  if( command == DACmd )
    { BeamOut->setDA(DACmd->GetNewIntValue(newValue));}
  if( command == DZCmd )
    { BeamOut->setDZ(DZCmd->GetNewIntValue(newValue));}
  if( command == ExCmd )
    { BeamOut->setEx(ExCmd->GetNewDoubleValue(newValue));}
  if( command == SrcCmd )
    { BeamOut->SetSource();}
  if( command == LvlCmd )
    { BeamOut->setLvlSchemeFile(newValue);}
  if( command == TACmd )
    { BeamOut->setTarA(TACmd->GetNewIntValue(newValue));}
  if( command == TZCmd )
    { BeamOut->setTarZ(TACmd->GetNewIntValue(newValue));}
  if( command == TExCmd )
    { BeamOut->setTarEx(TExCmd->GetNewDoubleValue(newValue));}
  if( command == TExFCmd )
    { BeamOut->setTFrac(TExFCmd->GetNewDoubleValue(newValue));}
  if( command == tauCmd )
    { BeamOut->settau(tauCmd->GetNewDoubleValue(newValue));}
/*  if( command == betaCmd )
    { BeamOut->setbeta(betaCmd->GetNewDoubleValue(newValue));}
  if( command == DoppZCmd )
    { BeamOut->setDoppZ(DoppZCmd->GetNewDoubleValue(newValue));}
  if( command == DoppYCmd )
    { BeamOut->setDoppY(DoppYCmd->GetNewDoubleValue(newValue));}
  if( command == DoppXCmd )
    { BeamOut->setDoppX(DoppXCmd->GetNewDoubleValue(newValue));}
*/
  if( command == RepCmd )
    { BeamOut->Report();}
  if( command == ThMTCmd )
    { BeamOut->SetThetaMaxTarget(ThMTCmd->GetNewDoubleValue(newValue));}
  if( command == AphTCmd )
    { BeamOut->SetAlphaTarget(AphTCmd->GetNewDoubleValue(newValue));}
  if( command == DistSigACmd )
    { BeamOut->SetThetaSigmaA(DistSigACmd->GetNewDoubleValue(newValue));}
  if( command == DistSigBCmd )
    { BeamOut->SetThetaSigmaB(DistSigBCmd->GetNewDoubleValue(newValue));}
  if( command == NQCmd )
    { BeamOut->SetNQ(NQCmd->GetNewIntValue(newValue));}
  if( command == SQCmd )
    { BeamOut->SelectQ(SQCmd->GetNewIntValue(newValue));}
  if( command == SCCmd )
    { BeamOut->SetQCharge(SCCmd->GetNewIntValue(newValue));}
  if( command == QUFCmd )
    { BeamOut->SetQUnReactedFraction(QUFCmd->GetNewDoubleValue(newValue));}
  if( command == QRFCmd )
    { BeamOut->SetQReactedFraction(QRFCmd->GetNewDoubleValue(newValue));}
  if( command == QKECmd )
    { BeamOut->SetQKE(QKECmd->GetNewDoubleValue(newValue));}
  if( command == QKEuCmd )
    { BeamOut->SetQKEu(QKEuCmd->GetNewDoubleValue(newValue));}
  if( command == a0cmd )
    { BeamOut->SetCoeff(0,a0cmd->GetNewDoubleValue(newValue));}
  if( command == a2cmd )
    { BeamOut->SetCoeff(2,a2cmd->GetNewDoubleValue(newValue));}
  if( command == a4cmd )
    { BeamOut->SetCoeff(4,a4cmd->GetNewDoubleValue(newValue));}
  if( command == a0Targetcmd )
    { BeamOut->SetTargetCoeff(0,a0Targetcmd->GetNewDoubleValue(newValue));}
  if( command == a2Targetcmd )
    { BeamOut->SetTargetCoeff(2,a2Targetcmd->GetNewDoubleValue(newValue));}
  if( command == a4Targetcmd )
    { BeamOut->SetTargetCoeff(4,a4Targetcmd->GetNewDoubleValue(newValue));}
}
