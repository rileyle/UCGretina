#include "PrimaryGeneratorAction_Messenger.hh"


PrimaryGeneratorAction_Messenger::PrimaryGeneratorAction_Messenger(PrimaryGeneratorAction* pga):PGA(pga) 
{ 
 
 
  PGADir = new G4UIdirectory("/Experiment/");
  PGADir->SetGuidance("Select the experimental conditions for simulations.");

  SrcDir = new G4UIdirectory("/Experiment/Source/");
  SrcDir->SetGuidance("Select source parameters for simulation.");

  RctDir = new G4UIdirectory("/Experiment/Reaction/");
  RctDir->SetGuidance("Select reaction parameters for simulation.");
  
  PGASCmd = new G4UIcmdWithoutParameter("/Experiment/RunSource",this);
  PGASCmd->SetGuidance("Select source calibration.");

  PGABCmd = new G4UIcmdWithoutParameter("/Experiment/RunBeam",this);
  PGABCmd->SetGuidance("Select in-beam simulations.");

  SrcCmd  = new G4UIcmdWithAString("/Experiment/Source/Set",this);
  SrcCmd->SetGuidance("Set source type (eu152, eu152_peaks (9 peaks, equal intensities), cs137, co56, co56_peaks (11 peaks, equal intensities), co60, ra226, ra226_peaks (10 peaks, equal intensities), photopeaks, au, white, background, bgwhite, simple, or muon)");
  SrcCmd->SetParameterName("Source type",false);
  SrcCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcECmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setEnergy",this);
  SrcECmd->SetGuidance("Set gamma-ray energy for the simple source type");
  SrcECmd->SetParameterName("Source Energy",false);
  SrcECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcWLECmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setWhiteLowE",this);
  SrcWLECmd->SetGuidance("Set lower gamma-ray energy limit for the white source");
  SrcWLECmd->SetParameterName("White Source Lower Energy Limit",false);
  SrcWLECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcWHECmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setWhiteHighE",this);
  SrcWHECmd->SetGuidance("Set upper gamma-ray energy limit for the white source");
  SrcWHECmd->SetParameterName("White Source Upper Energy Limit",false);
  SrcWHECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcXCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setX",this);
  SrcXCmd->SetGuidance("Set X position for the source");
  SrcXCmd->SetParameterName("Source X position",false);
  SrcXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcYCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setY",this);
  SrcYCmd->SetGuidance("Set Y position for the source");
  SrcYCmd->SetParameterName("Source Y position",false);
  SrcYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcZCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setZ",this);
  SrcZCmd->SetGuidance("Set Z position for the source");
  SrcZCmd->SetParameterName("Source Z position",false);
  SrcZCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcTFCmd = new G4UIcmdWithoutParameter("/Experiment/Source/OnTargetFace",this);
  SrcTFCmd->SetGuidance("Set source position on target face");
  
  SrcTBCmd = new G4UIcmdWithoutParameter("/Experiment/Source/OnTargetBack",this);
  SrcTBCmd->SetGuidance("Set source position on target back");

  SrcCollAngCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/CollimationAngle",this);
  SrcCollAngCmd->SetGuidance("Set angle of collimation (about Z) for the source");
  SrcCollAngCmd->SetParameterName("Collimation angle",false);
  SrcCollAngCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcCollDirCmd = new G4UIcmdWith3Vector("/Experiment/Source/CollimationDirection",this);
  SrcCollDirCmd->SetGuidance("Set direction of collimation for the source");
  SrcCollDirCmd->SetParameterName("X","Y","Z",false,false);
  SrcCollDirCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcRCmd = new G4UIcmdWithoutParameter("/Experiment/Source/Report",this);
  SrcRCmd->SetGuidance("Report source parameters");

  ROfCmd = new G4UIcmdWithoutParameter("/Experiment/Reaction/Off",this);
  ROfCmd->SetGuidance("Simulate only unreacted ions.");

  ROnCmd = new G4UIcmdWithoutParameter("/Experiment/Reaction/On",this);
  ROnCmd->SetGuidance("Simulate only Coulomb scattered ions.");

  SFrCmd = new G4UIcmdWithADouble("/Experiment/Reaction/UnscatteredFraction",this);
  SFrCmd->SetGuidance("Set fraction of Coulomb unscattered ions in the beam ");
  SFrCmd->SetParameterName("fraction",false);
  SFrCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
}



PrimaryGeneratorAction_Messenger::~PrimaryGeneratorAction_Messenger()
{
  
  delete PGADir;
  delete SrcDir;
  delete RctDir;
  delete PGASCmd;
  delete PGABCmd;
  delete SrcCmd;
  delete SrcECmd;
  delete SrcWLECmd;
  delete SrcWHECmd;
  delete SrcXCmd;
  delete SrcYCmd;
  delete SrcZCmd;
  delete SrcTFCmd;
  delete SrcTBCmd;
  delete SrcCollAngCmd;
  delete SrcCollDirCmd;
  delete SrcRCmd;
  delete ROnCmd;
  delete ROfCmd;
  delete SFrCmd;
}



void PrimaryGeneratorAction_Messenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 


  if( command == PGASCmd )
    {PGA ->SetSource();}

  if( command == PGABCmd )
    {PGA ->SetInBeam();}

  if( command == SrcCmd )
    {PGA ->SetSourceType(newValue);}

  if( command == SrcECmd )
    {PGA ->SetSourceEnergy(SrcECmd->GetNewDoubleValue(newValue));}

  if( command == SrcWLECmd )
    {PGA ->SetWhiteSourceLowE(SrcECmd->GetNewDoubleValue(newValue));}

  if( command == SrcWHECmd )
    {PGA ->SetWhiteSourceHighE(SrcECmd->GetNewDoubleValue(newValue));}

  if( command == SrcXCmd )
    {PGA ->SetSourceX(SrcXCmd->GetNewDoubleValue(newValue));}

  if( command == SrcYCmd )
    {PGA ->SetSourceY(SrcYCmd->GetNewDoubleValue(newValue));}

  if( command == SrcZCmd )
    {PGA ->SetSourceZ(SrcZCmd->GetNewDoubleValue(newValue));}

  if( command == SrcTFCmd )
    {PGA ->SetSourceOnTargetFace();}

  if( command == SrcTBCmd )
    {PGA ->SetSourceOnTargetBack();}

  if( command == SrcCollAngCmd )
    {PGA ->SetSourceCollAngle(SrcCollAngCmd->GetNewDoubleValue(newValue));}

  if( command == SrcCollDirCmd )
    {PGA ->SetSourceCollDirection(SrcCollDirCmd->GetNew3VectorValue(newValue));}

  if( command == SrcRCmd )
    {PGA ->SourceReport();}

  if( command == ROnCmd )
    {PGA ->ReactionOn();}

 if( command == ROfCmd )
    {PGA ->ReactionOff();}

 if( command == SFrCmd )
    {PGA ->SetFraction(SFrCmd->GetNewDoubleValue(newValue));}






}

