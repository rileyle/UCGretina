#ifndef CConvexPolyhedron_h
#define CConvexPolyhedron_h

//////////////////////////////////////////////////////////////////////////////////////////////
/// CConvexPolyhedron (D.B. 30/12/02)
///
/// Class description:
///
///  CConvexPolyhedron is a general Convex Polyhedron enclosed by flat surfaces.
///	It is built as a CSG solid
///	fPoints are the nPoints vertices used to define the planes
///	fPlanes the nPlanes enclosing planes
///	The faces can have any (>2) number of sides.
///	The coefficients of plane equations are obtained as an average
///	of the outside normals at each of the defining points.
///	This should avoid (or reduce??) problems of complanarity
///	originated by the precision of the fPoints
///
///  The class has been written (copied!) following the model of G4Trap (from Geant5.0)
///  G4Trap implements a 6-face convex polyhedron but treats the 2 planes
///  perpendicular to the z axis as a special case (which only complicates the code)
///
///  HepPolyhedron::createPolyhedron() method is used for visualization 
///  Faces with more than 4 sides are decomposed into 4-sided facets
///
///  To do:
///	   ~CConvexPolyhedron() to clear the allocated structures.
///	   More general constructors. The (real) one implemented so far assumes that the points
///	   are arranged in two surfaces joined by 4-sided faces (as normally used in AGATA).
///	   Accessors, Modifiers, ...
//////////////////////////////////////////////////////////////////////////////////////////////
 
#include "G4CSGSolid.hh"
#include "G4Point3D.hh"
#include "G4Plane3D.hh"
#include "G4Normal3D.hh"

#ifndef G4V10
#include "G4Point3DVector.hh"
#include "G4BoundingBox3D.hh"
//#else
//#include "G4Point3DList.hh"
#endif

class CConvexPolyhedron : public G4CSGSolid
{

  public:
    /////////////////////////////////////////////////////////////////
    /// copy constructor
    /////////////////////////////////////////////////////////////////
    CConvexPolyhedron(const CConvexPolyhedron& orig);

    /////////////////////////////////////////////////////////////////
    /// Constructor for MarsView: 2 faces of n/2 points connected by quadrangles
    /////////////////////////////////////////////////////////////////
#ifndef G4V10
    CConvexPolyhedron(const G4String& pName, const G4Point3DVector& pVec);
#else
    CConvexPolyhedron(const G4String& pName, const std::vector<G4Point3D>& pVec);
#endif
    /////////////////////////////////////////////////////////////////
    /// Generic Constructor: points and descriptors of faces
    /////////////////////////////////////////////////////////////////
#ifndef G4V10
    CConvexPolyhedron(const G4String& pName, const G4Point3DVector& pVec, 
                             const G4int nFaces, const std::vector<G4int>& theFaces);
#else
    CConvexPolyhedron(const G4String& pName, const std::vector<G4Point3D>& pVec, 
                             const G4int nFaces, const std::vector<G4int>& theFaces);
#endif
 
    ///////////////////////////////////////////////
    /// Destructor
    //////////////////////////////////////////////
    virtual ~CConvexPolyhedron() ;


  private:
    G4int       nPoints;            //> number of points
    G4Point3D  *fPoints;            //> the points as G4Point3D
    G4Point3D  *cCenter;            //> the center of the convex polyhedron (average of fPoints)

#ifndef G4V10
  private:
    G4BoundingBox3D* bbox;          //> used for a fast test in Inside()
#endif
    
  private:
    G4int        nPlanes;           //> number of enclosing planes
    G4int      **iPlanes;           //> descriptors of the planes (number of points, index of each point)
    G4Plane3D   *fPlanes;           //> the equation of the planes

  private:
    typedef double P3D[3];
    typedef int    FAC[4];
    P3D           *xyz;             //> the nPoints vertices for HepPolyhedron::createPolyhedron()
    int            nFacets;
    FAC           *fff;             //> the nFacets facets   for HepPolyhedron::createPolyhedron()

  public:
    inline G4int     GetnPoints ()         const {return nPoints;     };
    inline G4int     GetnPlanes ()         const {return nPlanes;     };
    inline G4Point3D GetPoints  (G4int nn)       {return fPoints[nn]; };
    inline G4Point3D* GetCenter ()               { return cCenter;    };
    

  public:
    G4bool MovePlane(const G4int nn, const G4double dist);   //> shifts the plane along its normal (+ --> outwards)
    G4bool CalculateExtent( const EAxis pAxis, const G4VoxelLimits& pVoxelLimit,
                                  const G4AffineTransform& pTransform, G4double& pMin, G4double& pMax ) const;    
        
  public:
    EInside Inside( const G4ThreeVector& p ) const;
         
  public:
    G4ThreeVector SurfaceNormal( const G4ThreeVector& p ) const;

  public:
    G4double DistanceToIn( const G4ThreeVector& p, const G4ThreeVector& v) const;
    G4double DistanceToIn( const G4ThreeVector& p ) const;

  public:
    G4double DistanceToOut( const G4ThreeVector& p ) const;
    G4double DistanceToOut( const G4ThreeVector& p, const G4ThreeVector& v,
                                  const G4bool calcNorm=false, G4bool *validNorm=0, G4ThreeVector *n=0) const;

  public:
    G4GeometryType GetEntityType() const;
    std::ostream& StreamInfo( std::ostream& os ) const;

  //////////////////////////////////////////////////////////////////
  /// Visualisation functions
  ///////////////////////////////////////////////////////////////////
  public:
    void          DescribeYourselfTo ( G4VGraphicsScene& scene  ) const;
    G4Polyhedron* CreatePolyhedron   () const;
#ifndef G4V10
    G4NURBS*      CreateNURBS        () const;
#endif

  protected:
    G4ThreeVectorList* CreateRotatedVertices( const G4AffineTransform& pTransform ) const ;

  private:
    void      MakePlane(G4int ii);  //> construct equation of plane from its defining points
    G4bool    MakePoint(G4int ii);  //> determine point as intersection of 3 planes
    void      CheckConvexity(const G4double dmax);
    void      PreparePolyhedron();
    G4Point3D SolvePointEq(G4Plane3D& pl0, G4Plane3D& pl1, G4Plane3D& pl2);

};

#endif
