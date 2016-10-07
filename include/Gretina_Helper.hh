#ifndef AgataAncillaryHelper_h
#define AgataAncillaryHelper_h 1


#include "globals.hh"
#include "G4Point3D.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include <vector>

using namespace std;

class CConvexPolyhedron;
class G4Tubs;
class G4Polycone;
class G4IntersectionSolid;
class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;

////////////////////////////////////////////////////////////////////////////////////
/// This class handles the information about the encapsulated detectors of AGATA
/// (size, ...) together with the objects actually generating the detectors
/// within the simulation. It is basically a struct, but was written as a class
/// where all members are public since the compiler does not accept to allocate
/// a std::vector with a structure. It is also used for the cryostats.
////////////////////////////////////////////////////////////////////////////////////
class CpolyhPoints
{
  // default empty creator/destructor
  public:
     CpolyhPoints() {};
    ~CpolyhPoints() {};
  
  ///////////////////////////////
  /// All members are public!!!  
  ///////////////////////////////
  public:  
    G4int  whichGe;        //> kind of crystal (for the cryostats: cluster they belong to)
    G4int  whichCrystal;   //> for the cryostats: crystal they belong to
    G4int  whichWall;      //> for the cryostats: kind of part composing the cryostat
    
  public:  
    G4int  npoints;        //> number of points composing the solid
    G4int  nfaces;         //> number of faces composing the solid
    
  public:  
    std::vector<G4Point3D> vertex;  //> 3D-vertexes of the solid
    std::vector<G4int>     ifaces;  //> points (of vertex) forming each face
    
  ////////////////////
  /// Cylinder data
  ///////////////////
  public:  
    G4double  tubr, tubR, tubL;      //> inner and outer radius, length
    G4double  tubX, tubY, tubZ;      //> position of the cylinder with respect to the polyhedron (for the intersection)
    G4double  thick;                 //> distance from the front face to the coaxial hole
    
  public:
    G4bool    cylinderMakesSense;    //> false: dimensions are not compatible with a real coaxial crystal
                                     //>       (thus the cylinder will not be built)  

  ///////////////////
  /// Capsule data
  //////////////////
  public:
    G4double  capThick;              //> thickness
    G4double  capSpace;              //> crystal-capsule spacing

  public:
    G4bool    makeCapsule;           //> true: encapsulation will be built

  ///////////////////////////
  /// passivated areas size
  ///////////////////////////
  public:
    G4double  passThick1;            //> at the back of the crystal
    G4double  passThick2;            //> around the coaxial hole
    
  public:
    G4Point3D centerFace1;           //> center of front face
    G4Point3D centerFace2;           //> center of back face  
    
  public:  
    G4double  zFace1;                //> z-coordinate of the center of the front face
    G4double  zFace2;                //> z-coordinate of the center of the back face
    G4double  zCenter;               //> z-coordinate of the center of the crystal
    
  ///////////////////
  /// segmentation
  ///////////////////
  public:  
    G4int                  nslice;   //> number of slices in which the crystal is divided
    std::vector<G4double>  zSliceI;  //> z-coordinate of the slices at the crystal axis
    std::vector<G4double>  zSliceO;  //> z-coordinate of the slices at the outer surface
    
  public:  
    G4double  minR;                  //> minimum radius of a cylinder surrounding the polyhedron
    
  public:  
    G4double  colx, coly, colz;      //> RGB components defining the colour

  ///////////////////
  /// planars
  ///////////////////
  public: 
    G4double  guardThick[4];         //> thickness of the guardring
    G4bool    isPlanar;              //> true: treat detector as a planar 

  public:
    G4double  segSize_x;             //> x-size of the segments
    G4double  segSize_y;             //> y-size of the segments
    G4double  maxSize_x;             //> max x-coordinate of the crystal
    G4double  maxSize_y;             //> max y-coordinate of the crystal
    G4double  minSize_x;             //> min x-coordinate of the crystal
    G4double  minSize_y;             //> min y-coordinate of the crystal
    G4int     nSeg_x;                //> number of pixels in x direction
    G4int     nSeg_y;                //> number of pixels in y direction
     
    
  /////////////////////////////////////////////////////
  /// Objects defining the detector in the simulation
  /////////////////////////////////////////////////////
  public:  
    CConvexPolyhedron   *pPoly;      //> original polyhedron
    G4Polycone          *pCoax;      //> cylinder 
    G4IntersectionSolid *pCaps;      //> their intersection
    G4LogicalVolume     *pDetL;      //> its logical
    
  public:  
    G4Polycone          *pTubs1;     //> cylinder (passive area behind detector)     
    G4IntersectionSolid *pCaps1;     //> intersection with polyhedron		    
    G4LogicalVolume     *pDetL1;     //> its logical				    
    G4VPhysicalVolume   *pDetP1;     //> passivated area (back)			    
   
  public:  
    G4Polycone          *pCoax2;     //> cylinder (passive area at the coaxial hole) 
    G4IntersectionSolid *pCaps2;     //> intersection with polyhedron		    
    G4LogicalVolume     *pDetL2;     //> its logical				    
    G4VPhysicalVolume   *pDetP2;     //> passivated area (coax)			    
    
  public:  
    G4VisAttributes     *pDetVA;     //> visualization attributes
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This "structure" (written as a class as explained above) handles the transformations needed to
/// place the AGATA clusters into space.
////////////////////////////////////////////////////////////////////////////////////////////////////
class CeulerAngles
{
  // default empty creator/destructor
  public:
     CeulerAngles() {};
    ~CeulerAngles() {};

  // all members are public  
  public:  
    G4int     whichGe;       //> type of cluster which should be placed   
    G4int     numPhys;       //> number of cluster which should be placed

  public:  
    G4double  ps, th, ph;    //> Rotation angles; they are needed since one cannot easily extract
                             //> them from the rotation matrix

  public:
    G4ThreeVector trasl;     //> traslation

  public:
    G4RotationMatrix rotMat; //> rotation matrix

  public:
    G4Transform3D *pTransf;  //> full transformation
};

class G4AssemblyVolume;
////////////////////////////////////////////////////////////////////////////////////////////////////
/// This "structure" (written as a class as explained above) handles the transformations needed to
/// place the crystals within an AGATA cluster
////////////////////////////////////////////////////////////////////////////////////////////////////
class CclusterAngles
{
  // default empty creator/destructor
  public:
     CclusterAngles() {};
    ~CclusterAngles() {};

  // all members are public  
  public:
    G4int                      whichClus;   //> type of cluster  
    
  public:
    G4int                      nsolids;     //> number of crystals arranged to form a cluster 
    std::vector<CeulerAngles>  solids;      //> each component of this array contains the information
                                            //> needed to place a single crystal within the cluster
    
  public:
    G4int                      nwalls;      //> number of parts composing the cryostat  
#ifdef ANTIC
    G4int                      nantic;      //> number of parts composing the antiCompton 
#endif
    
  public:
    G4AssemblyVolume*          pAssV;       //> pointer to the G4AssemblyVolume forming the cluster
};

#endif
