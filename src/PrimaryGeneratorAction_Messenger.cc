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
  SrcWLECmd->SetGuidance("Set lower gamma-ray energy limit for the white/bgwhite source");
  SrcWLECmd->SetParameterName("White Source Lower Energy Limit",false);
  SrcWLECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcWHECmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setWhiteHighE",this);
  SrcWHECmd->SetGuidance("Set upper gamma-ray energy limit for the white/bgwhite source");
  SrcWHECmd->SetParameterName("White Source Upper Energy Limit",false);
  SrcWHECmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcMultCmd = new G4UIcmdWithAnInteger("/Experiment/Source/setMultiplicity",this);
  SrcMultCmd->SetGuidance("Set multiplicity for the white/bgwhite source");
  SrcMultCmd->SetParameterName("White Source Multiplicity",false);
  SrcMultCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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

  SrcRCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setR",this);
  SrcRCmd->SetGuidance("Set the radius of the source disk");
  SrcRCmd->SetParameterName("Source radius",false);
  SrcRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcDXCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setDX",this);
  SrcDXCmd->SetGuidance("Set the width of the source in the  nondispersive direction");
  SrcDXCmd->SetParameterName("Source width",false);
  SrcDXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcDYCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setDY",this);
  SrcDYCmd->SetGuidance("Set the height of the source in the dispersive direction");
  SrcDYCmd->SetParameterName("Source height",false);
  SrcDYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcSigXCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setSigmaX",this);
  SrcSigXCmd->SetGuidance("Set the width of the Gaussian source distribution in the nondispersive direction");
  SrcSigXCmd->SetParameterName("Source Gaussian width",false);
  SrcSigXCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcSigYCmd = new G4UIcmdWithADoubleAndUnit("/Experiment/Source/setSigmaY",this);
  SrcSigYCmd->SetGuidance("Set the height of the Gaussian source distribution in the dispersive direction");
  SrcSigYCmd->SetParameterName("Source Gaussian height",false);
  SrcSigYCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

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

  SrcThetaFileCmd = new G4UIcmdWithAString("/Experiment/Source/ThetaFile",this);
  SrcThetaFileCmd->SetGuidance("Set the name of the file specifying the theta distribution of the emitted particles.");
  SrcThetaFileCmd->SetParameterName("choice",false);
  SrcThetaFileCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcRepCmd = new G4UIcmdWithoutParameter("/Experiment/Source/Report",this);
  SrcRepCmd->SetGuidance("Report source parameters");

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
  delete SrcMultCmd;
  delete SrcXCmd;
  delete SrcYCmd;
  delete SrcZCmd;
  delete SrcRCmd;
  delete SrcDXCmd;
  delete SrcDYCmd;
  delete SrcSigXCmd;
  delete SrcSigYCmd;
  delete SrcTFCmd;
  delete SrcTBCmd;
  delete SrcCollAngCmd;
  delete SrcCollDirCmd;
  delete SrcThetaFileCmd;
  delete SrcRepCmd;
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

  if( command == SrcMultCmd )
    {PGA ->SetWhiteSourceMult(SrcMultCmd->GetNewIntValue(newValue));}

  if( command == SrcXCmd )
    {PGA ->SetSourceX(SrcXCmd->GetNewDoubleValue(newValue));}

  if( command == SrcYCmd )
    {PGA ->SetSourceY(SrcYCmd->GetNewDoubleValue(newValue));}

  if( command == SrcZCmd )
    {PGA ->SetSourceZ(SrcZCmd->GetNewDoubleValue(newValue));}

  if( command == SrcRCmd )
    {PGA ->SetSourceR(SrcRCmd->GetNewDoubleValue(newValue));}

  if( command == SrcDXCmd )
    {PGA ->SetSourceDX(SrcDXCmd->GetNewDoubleValue(newValue));}

  if( command == SrcDYCmd )
    {PGA ->SetSourceDY(SrcDYCmd->GetNewDoubleValue(newValue));}

  if( command == SrcSigXCmd )
    {PGA ->SetSourceSigmaX(SrcSigXCmd->GetNewDoubleValue(newValue));}

  if( command == SrcSigYCmd )
    {PGA ->SetSourceSigmaY(SrcSigYCmd->GetNewDoubleValue(newValue));}

  if( command == SrcTFCmd )
    {PGA ->SetSourceOnTargetFace();}

  if( command == SrcTBCmd )
    {PGA ->SetSourceOnTargetBack();}

  if( command == SrcCollAngCmd )
    {PGA ->SetSourceCollAngle(SrcCollAngCmd->GetNewDoubleValue(newValue));}

  if( command == SrcCollDirCmd )
    {PGA ->SetSourceCollDirection(SrcCollDirCmd->GetNew3VectorValue(newValue));}

  if( command == SrcThetaFileCmd )
    {PGA -> SetSourceThetaFile(newValue);}
  
  if( command == SrcRepCmd )
    {PGA ->SourceReport();}

}

