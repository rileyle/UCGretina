#ifndef Gretina_Array_h
#define Gretina_Array_h 1

#include "G4VUserDetectorConstruction.hh"
#include "DetectorConstruction.hh"
#include "TrackerGammaSD.hh"
#include "RunAction.hh"
#include "CConvexPolyhedron.hh"
#include "Gretina_Helper.hh"

#include "G4AssemblyVolume.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"

#include "globals.hh"
#include "G4Point3D.hh"
#include "G4Plane3D.hh"
#include "G4Normal3D.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include <vector>

#include "G4SystemOfUnits.hh"

class Gretina_Array_Messenger;

using namespace std;
using namespace CLHEP;

/////////////////////////////////////////////////////////////////////
/// This class handles the construction of the actual GRETINA array
////////////////////////////////////////////////////////////////////
class Gretina_Array
{
  public:
    Gretina_Array();
    ~Gretina_Array();
    
  private:
    Gretina_Array_Messenger *myMessenger;
    
  /////////////////////////////////////////////////
  /// Files from which actual geometry is read 
  ///////////////////////////////////////////////// 
  private:
    G4String                    iniPath;     //> directory where the files are located
    G4String                    eulerFile;   //> angles and positions to place the clusters into space
    G4String                    solidFile;   //> shape of the crystals
    G4String                    sliceFile;   //> segmentation
    G4String                    wallsFile;   //> cryostats
    G4String                    clustFile;   //> arrangement of the crystals and cryostats within a cluster
    
  private:
    G4String                    directoryName;  //> for the command line  
    
  /////////////////////////////////////
  /// materials (pointers and names)
  ///////////////////////////////////// 
  private:
    G4Material                  *matCryst;       //> crystals
    G4Material                  *matWalls;       //> encapsulation
    G4Material                  *matBackWalls;   //> behind the crystals
    G4Material                  *matHole;        //> vacuum within the cryostat
    G4Material                  *matCryo;        //> cryostats

  private:
    G4String                    matCrystName;     //> crystals
    G4String                    matWallsName;     //> cryostats and encapsulation
    G4String                    matBackWallsName; //> behind the crystals
    G4String                    matHoleName;      //> vacuum within the cryostat
    G4String                    matCryoName;      //> cryostats
     
  ///////////////////////////////////////////////////////////////////////
  /// structures needed to store geometry data during the construction  
  ///////////////////////////////////////////////////////////////////////
  private:
    std::vector<CeulerAngles>   euler;        //> angles and positions to place the clusters into space
    std::vector<CpolyhPoints>   pgons;        //> shape of the crystals
    std::vector<CpolyhPoints>   walls;        //> cryostats
    std::vector<CclusterAngles> clust;        //> arrangement of the crystals and cryostats within a cluster
    std::vector<CpolyhPoints>   capsO;        //> encapsulation (outer size)
    std::vector<CpolyhPoints>   capsI;        //> encapsulation (inner size)

  private:
    G4int                       nEuler;       //> number of clusters composing the array
    G4int                       nPgons;       //> number of different crystal shapes within the array
    G4int                       nClAng;       //> number of crystals composing a cluster
    G4int                       nWalls;       //> number of cryostat parts within a cluster
    G4int                       maxPgons;     //> maximum index of crystal shapes
    G4int                       nDets;
    G4int                       nClus;
    G4int                       iCMin;
    G4int                       iCMax;
    G4int                       iGMin;
    G4int                       iGMax;
    G4int                       maxSec;
    G4int                       maxSli;

  private:
    std::vector<G4int>          crystType;    //> lookup table detector number --> crystal shape
    std::vector<G4int>          planarLUT;    //> lookup table detector number --> planar or not

  private:
    G4int                       nWlTot;       //> total number of cryostat parts within the array
    G4int                       maxSolids;    //> maximum number of solids within a cluster
    
  /////////////////////////////////////////////////
  /// structures needed to build the segmentation
  ////////////////////////////////////////////////
  private:
    std::vector<CpolyhPoints>   pgSegLl;      //> segments on lower Left  side of edges
    std::vector<CpolyhPoints>   pgSegLu;      //> segments on upper Left  side of edges
    std::vector<CpolyhPoints>   pgSegRl;      //> segments on lower Right side of edges
    std::vector<CpolyhPoints>   pgSegRu;      //> segments on upper Right side of edges

  private:
    std::vector<G4int>          nSegments;
    std::vector<G4int>          tSegments;
    G4int                       totSegments;
    G4int                       nSeg;
    
  private:
    std::vector<G4double>       segVolume;    //> the volume of the (composite) segment
    std::vector<G4Point3D>      segCenter;    //> the center of mass of the (composite) segment

  private:
    G4int                       stepFactor;    //> integration step for the calculation of segment volume
    G4bool                      stepHasChanged;//> true: integration step was changed and segment volumes
                                               //>       should be recomputed
   
  ////////////////////////////////////////////
  /// size of the equivalent germanium shell
  ////////////////////////////////////////////
  private:
    G4double                    arrayRmin;     //> inner radius
    G4double                    arrayRmax;     //> outer radius
    
  ///////////////////////////////////////////
  /// rotation applied to the whole array
  //////////////////////////////////////////
  private:
    G4double                    thetaShift;   //> theta
    G4double                    phiShift;     //> phi
    G4double                    thetaPrisma;  //> thetaPrisma

  ///////////////////////////////////////////
  /// traslation applied to the whole array
  //////////////////////////////////////////
  private:
    G4ThreeVector               posShift;   

  ///////////////////////////////////////////                               //LR
  /// Cryostats             
  //////////////////////////////////////////
  public:
    void SetCryostats(){cryostatStatus = true;};                            //LR

  private:
    G4ThreeVector               cryostatPos0;                               //LR
    G4ThreeVector               cryostatPos;                                //LR
    G4RotationMatrix            cryostatRot;                                //LR
  //G4double                    cryostatRadius;                             //LR
  //G4double                    cryostatLength;                             //LR
    G4double                    cryostatZplanes[7];                         //LR
    G4double                    cryostatRinner[7];                          //LR
    G4double                    cryostatRouter[7];                          //LR
  /////////////////
  /// some flags 
  //////////////// 
  private:
    G4bool                      usePassive;     //> true: passive areas of the crystals will be generated
    G4bool                      drawReadOut;    //> true: segments will be visualized
    G4bool                      useAncillary;   //> true: ancillary detectors will be constructed
    G4bool                      makeCapsule;    //> true: encapsulation will be generated
    G4bool                      useCylinder;    //> true: the intersection with the cylinder will be considered
    G4bool                      cryostatStatus; //> true: include cryostats behind clusters //LR
    G4bool                      readOut;        //> true: a segmentation have been defined  
    G4bool                      printVolumes;   //> true: report crystal, segment, and passive volumes
  
  //////////////////////////////
  ///////// Methods  ///////////
  //////////////////////////////
  private:
    void     InitData();

  //////////////////////////
  /// read the input files
  /////////////////////////
  private:
    void      ReadEulerFile();
    void      ReadSolidFile();
    void      ReadSliceFile();
    void      ReadWallsFile();
    void      ReadClustFile();
  
  //////////////////////////////////////////////////////////////
  /// look for the materials starting from the material names
  /////////////////////////////////////////////////////////////
  private:
    G4int     FindMaterials();
    
  /////////////////////////////////////////////////////////
  /// Construct the various elements composing the array
  /////////////////////////////////////////////////////////  
  private:
    void      ConstructGeCrystals  ();    
    void      ConstructTheCapsules ();    
    void      ConstructTheClusters ();    
    void      ConstructTheWalls    ();
    
  ////////////////////////////////
  /// placement of the elements  
  ////////////////////////////////
  private:
    void      PlaceTheClusters     ();

  //////////////////////////////////
  /// Construction of the segments
  //////////////////////////////////
  
  private:
    //////////////////////////////
    /// Construct the segments
    //////////////////////////////
    void      ConstructSegments        ();
    /////////////////////////////////////////////////////////////////////////////////
    /// Calculate the vertexes of the segments starting from the original polyhedra
    /////////////////////////////////////////////////////////////////////////////////
    G4int     CalculateSegments        ( G4int );
    ///////////////////////////////////
    /// Checks for possible overlaps
    //////////////////////////////////
    G4int     CheckOverlap             ( G4int, G4int, G4int );
    /////////////////////////////////////////////////////////////////////////////////////////
    /// Calculates volume and center of the segments (each of them composed of more parts!)
    //////////////////////////////////////////////////////////////////////////////////////////
    void      CalculateVolumeAndCenter ( G4int, G4int, G4int, G4double );
    ////////////////////////////////////////////////////////////////////////////////////////
    /// Calculates the intersection between a plane and a line passing through two points
    ///////////////////////////////////////////////////////////////////////////////////////
    G4Point3D XPlaneLine               ( const G4Plane3D &vv, const G4Point3D &pA,  const G4Point3D &pB );
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Calculates segment number (corresponding to a given detector and position relative to the crystal)
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    G4int     GetCoaxSegmentNumber         ( G4int, G4ThreeVector );      

  //////////////////////////////////
  /// writes out the information
  /////////////////////////////////
  private:
    void WritePositions               ( std::ofstream &outFileLMD, G4double=1.*mm );
    void WriteSegmentPositions        ( std::ofstream &outFileLMD, G4double=1.*mm );
    void WriteCrystalPositions        ( std::ofstream &outFileLMD, G4double=1.*mm );
    void WriteCrystalTransformations  ( std::ofstream &outFileLMD, G4double=1.*mm );
  
  ////////////////
  /// placement
  ///////////////
  public:
    void Placement ();

  //////////////////////////////////////////////
  /// public interface to the private method!
  ////////////////////////////////////////////////
  public:
    G4int     GetSegmentNumber         ( G4int, G4ThreeVector );    
  
  //////////////////////////////////
  /// writes out the information
  /////////////////////////////////
  public:
    void WriteHeader                  ( std::ofstream &outFileLMD, G4double=1.*mm );
    void WriteSegmentAngles           ( G4String, G4int = 0 );
    void WriteCrystalAngles           ( G4String );

  //////////////////////////////////////////////////
  ///////////// public methods for the messenger
  ////////////////////////////////////////////////// 
  public:       
    void SetSolidFile           ( G4String );
    void SetWallsFile           ( G4String );
    void SetAngleFile           ( G4String );
    void SetSliceFile           ( G4String );
    void SetClustFile           ( G4String );

  public:       
    void SetDetMate             ( G4String );
    void SetWallsMate           ( G4String );
    void SetBackWallsMate       ( G4String );

  public:      
    void SetThetaShift          ( G4double );
    void SetPhiShift            ( G4double );
    void SetPosShift            ( G4ThreeVector );
    void SetThetaPrisma         ( G4double );

  public:      
    void SetMakeCapsules        ( G4bool );
    void SetUseCylinder         ( G4bool );
    void SetUsePassive          ( G4bool );
    void SetDrawReadOut         ( G4bool );
    void SetWriteSegments       ( G4bool );
    void SetPrintVolumes        ( G4bool flag){ printVolumes = flag; }

public:       
    void SetStep                ( G4int  ); 
     
  public:
    void ShowStatus ();
    
  ///////////////////////////////////////////
  //////////////// inline "get" methods
  ///////////////////////////////////////////
  public:
    inline G4double              GetThetaShift      () { return thetaShift;  };
    inline G4double              GetPhiShift        () { return phiShift;    };
    
  public:
    inline G4String              GetEulerFile       () { return eulerFile;   };
    inline G4String              GetSolidFile       () { return solidFile;   };
    inline G4String              GetSliceFile       () { return sliceFile;   };
    inline G4String              GetWallsFile       () { return wallsFile;   };
    
  public:
    inline G4bool                GetDrawReadOut     () { return drawReadOut; };
    inline G4bool                GetReadOut         () { return readOut; };
};

#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWith3Vector;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

class Gretina_Array_Messenger: public G4UImessenger
{
  public:
    Gretina_Array_Messenger(Gretina_Array*);
   ~Gretina_Array_Messenger();
    
  private:
    Gretina_Array*        myTarget;
    
  private:
    G4UIcmdWithAString*        SetSolidCmd;
    G4UIcmdWithAString*        SetAngleCmd;
    G4UIcmdWithAString*        SetWallsCmd;
    G4UIcmdWithAString*        SetClustCmd;
    G4UIcmdWithAString*        SetSliceCmd;
    G4UIcmdWithAString*        DetMatCmd;
    G4UIcmdWithAString*        WalMatCmd;
    G4UIcmdWithAString*        BackWalMatCmd;
    G4UIcmdWithAString*        RotateArrayCmd;
    G4UIcmdWithADouble*        RotatePrismaCmd;
    G4UIcmdWith3Vector*        TraslateArrayCmd;
    G4UIcmdWithAString*        WriteAnglesCmd;
    G4UIcmdWithAString*        WriteCryAnglesCmd;
    G4UIcmdWithABool*          EnableCylCmd;
    G4UIcmdWithoutParameter*   cryostatCmd;
    G4UIcmdWithABool*          DisableCylCmd;
    G4UIcmdWithABool*          EnablePassiveCmd;
    G4UIcmdWithABool*          DisablePassiveCmd;
    G4UIcmdWithABool*          DontDrawReadOutCmd;
    G4UIcmdWithABool*          DrawReadOutCmd;
    G4UIcmdWithABool*          EnableCapsulesCmd;
    G4UIcmdWithABool*          DisableCapsulesCmd;            
    G4UIcmdWithAnInteger*      SetStepCmd;
    G4UIcmdWithoutParameter*   printVolCmd;

public:
    void SetNewValue(G4UIcommand*, G4String);
};

#endif
