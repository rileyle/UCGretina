#include "CConvexPolyhedron.hh"
#include "G4ThreeVector.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include <iomanip>

using namespace std; // needed for swap()

#ifdef G4V10
using namespace CLHEP;
#endif


CConvexPolyhedron::~CConvexPolyhedron()
{
  if(fPoints) delete [] fPoints;
  if(cCenter) delete    cCenter;
#ifndef G4V10
  if(bbox)    delete    bbox;
#endif
  if(fPlanes) delete [] fPlanes;
  for(G4int nn = 0; nn < nPlanes; nn++ )
    if(iPlanes[nn]) delete [] iPlanes[nn];
  if(iPlanes) delete [] iPlanes;
  if(xyz)     delete [] xyz;
  if(fff)     delete [] fff;
}

// the copy constructor
#ifndef G4V10
CConvexPolyhedron::CConvexPolyhedron(const CConvexPolyhedron& orig)
  : G4CSGSolid(orig.GetName()), nPoints(0), fPoints(0), cCenter(0), bbox(0),
                                nPlanes(0), iPlanes(0), fPlanes(0),
                                xyz(0), nFacets(0), fff(0)
#else
CConvexPolyhedron::CConvexPolyhedron(const CConvexPolyhedron& orig)
  : G4CSGSolid(orig.GetName()), nPoints(0), fPoints(0), cCenter(0),
                                nPlanes(0), iPlanes(0), fPlanes(0),
                                xyz(0), nFacets(0), fff(0)
#endif
{
  nPoints = orig.nPoints;
  fPoints = new G4Point3D[nPoints];
  cCenter = new G4Point3D;
#ifndef G4V10
  bbox    = new G4BoundingBox3D;
#endif
  
  G4int nn;
  for( nn = 0; nn < nPoints; nn++ )
    fPoints[nn] = orig.fPoints[nn];

  cCenter = orig.cCenter;
#ifndef G4V10 
  bbox    = orig.bbox;
#endif

  nPlanes = orig.nPlanes;
  iPlanes = new G4int*[nPlanes];
  fPlanes = new G4Plane3D[nPlanes];
  
  for( nn = 0; nn < nPlanes; nn++ ) {
    G4int  nSides = *orig.iPlanes[nn];      // how many points
    iPlanes[nn] = new int[1+nSides];
    memcpy(iPlanes[nn], orig.iPlanes[nn], (nSides+1)*sizeof(G4int));
    fPlanes[nn] = orig.fPlanes[nn];
  }

  xyz = new P3D[nPoints];
  memcpy(xyz, orig.xyz, nPoints*sizeof(P3D));

  nFacets = orig.nFacets;
  fff = new FAC[nFacets];
  memcpy(fff, orig.fff, nFacets*sizeof(FAC));

}

// constructor for 2 polygonal faces of n/2 points with quadrangular sides
#ifndef G4V10
CConvexPolyhedron::CConvexPolyhedron(const G4String& pName, const G4Point3DVector &pVec)
  : G4CSGSolid(pName), nPoints(0), fPoints(0), cCenter(0), bbox(0),
                       nPlanes(0), iPlanes(0), fPlanes(0),
                       xyz(0), nFacets(0), fff(0)
#else
    CConvexPolyhedron::CConvexPolyhedron(const G4String& pName, const std::vector<G4Point3D> &pVec)
  : G4CSGSolid(pName), nPoints(0), fPoints(0), cCenter(0),
                       nPlanes(0), iPlanes(0), fPlanes(0),
                       xyz(0), nFacets(0), fff(0)
#endif
{
  nPoints = pVec.size();

  if(nPoints < 6 || nPoints%2) {
    G4cout << "ERROR - CConvexPolyhedron::CConvexPolyhedron(): " << GetName() << G4endl
           << "        Invalid number of vertices: " << nPoints << G4endl;
    G4cerr << "ERROR - CConvexPolyhedron::CConvexPolyhedron(): " << GetName() << G4endl
           << "        Invalid number of vertices: " << nPoints << G4endl;
#if defined G4V495 || G4V496 || G4V10
    G4Exception("CConvexPolyhedron::CConvexPolyhedron()","Error2",FatalException,"Invalid number of vertices");
#else
    G4Exception("CConvexPolyhedron::CConvexPolyhedron() - Invalid number of vertices");
#endif
  }

  G4int nSides = nPoints/2;
  fPoints = new G4Point3D[nPoints];
  cCenter = new G4Point3D;
#ifndef G4V10
  bbox    = new G4BoundingBox3D;
#endif

  G4int ii;
  for( ii = 0; ii < nPoints; ++ii ) {
    fPoints[ii] = pVec[ii];
    *cCenter += fPoints[ii];
  }
  *cCenter /= nPoints;

#ifndef G4V10
  bbox->Init(*cCenter);
  for( ii = 0; ii < nPoints; ++ii ) {
    bbox->Extend(fPoints[ii]);
  }
#endif
  
  nPlanes = nSides + 2;
  iPlanes = new G4int*[nPlanes];
  fPlanes = new G4Plane3D[nPlanes];

  iPlanes[0] = new int[1+nSides];
  iPlanes[0][0] = nSides;
  for ( ii = 0; ii < nSides; ++ii )
    iPlanes[0][ii+1] = ii;

  iPlanes[1] = new int[1+nSides];
  iPlanes[1][0] = nSides;
  for ( ii = 0; ii < nSides; ++ii )
    iPlanes[1][ii+1] = nSides + ii;

  for (ii = 0; ii < nSides; ++ii) {
    iPlanes[2+ii]    = new G4int[5];
    iPlanes[2+ii][0] = 4;
    iPlanes[2+ii][1] =           ii;
    iPlanes[2+ii][2] = nSides +  ii;
    iPlanes[2+ii][3] = nSides + (ii+1)%nSides;
    iPlanes[2+ii][4] =          (ii+1)%nSides;
  } 

  for( ii = 0; ii < nPlanes; ii++ )
    MakePlane(ii);

//CheckConvexity(kCarTolerance/2.);
// given that we read the points from a text file, this is a more realistic precision
  CheckConvexity(0.05*micrometer);

  G4Point3D cc(0., 0., 0.);
  for( ii = 0; ii < nPoints; ii++ ) {   // recalculate the points
    MakePoint(ii);
    cc += fPoints[ii];
#ifndef G4V10
    bbox->Extend(fPoints[ii]);
#endif
  }
  cc /= nPoints;
  *cCenter = cc;

#ifndef G4V10
  bbox->Init(cc);
  for( ii = 0; ii < nPoints; ii++ ) {   // adjust bbox
    bbox->Extend(fPoints[ii]);
  }
#endif
  
  CheckConvexity(kCarTolerance/2.);   // now we can check with the default precision

  PreparePolyhedron();

}

// generic constructor: the points and descriptors of the faces
#ifndef G4V10
CConvexPolyhedron::CConvexPolyhedron(const G4String& pName, const G4Point3DVector &pVec,
                                     const G4int nFaces, const std::vector<G4int> &theFaces)
  : G4CSGSolid(pName), nPoints(0), fPoints(0), cCenter(0), bbox(0),
                       nPlanes(0), iPlanes(0), fPlanes(0),
                       xyz(0), nFacets(0), fff(0)
#else
    CConvexPolyhedron::CConvexPolyhedron(const G4String& pName, const std::vector<G4Point3D> &pVec,
                                     const G4int nFaces, const std::vector<G4int> &theFaces)
  : G4CSGSolid(pName), nPoints(0), fPoints(0), cCenter(0),
                       nPlanes(0), iPlanes(0), fPlanes(0),
                       xyz(0), nFacets(0), fff(0)
#endif
{
  nPoints = pVec.size();

  if(nPoints < 4) {
    G4cout << "ERROR - CConvexPolyhedron::CConvexPolyhedron(): " << GetName() << G4endl
           << "        Invalid number of vertices: " << nPoints << G4endl;
    G4cerr << "ERROR - CConvexPolyhedron::CConvexPolyhedron(): " << GetName() << G4endl
           << "        Invalid number of vertices: " << nPoints << G4endl;
#if defined G4V495 || G4V496 || G4V10
    G4Exception("CConvexPolyhedron::CConvexPolyhedron()","Error3",FatalException,"Invalid number of vertices");
#else
    G4Exception("CConvexPolyhedron::CConvexPolyhedron() - Invalid number of vertices");
#endif
  }

  if(nFaces < 4) {
    G4cout << "ERROR - CConvexPolyhedron::CConvexPolyhedron(): " << GetName() << G4endl
           << "        Invalid number of faces: " << nFaces << G4endl;
    G4cerr << "ERROR - CConvexPolyhedron::CConvexPolyhedron(): " << GetName() << G4endl
           << "        Invalid number of faces: " << nFaces << G4endl;
#if defined G4V495 || G4V496 || G4V10
    G4Exception("CConvexPolyhedron::CConvexPolyhedron()","Error4",FatalException,"Invalid number of faces");
#else
    G4Exception("CConvexPolyhedron::CConvexPolyhedron() - Invalid number of faces");
#endif
  }

  fPoints = new G4Point3D[nPoints];
  cCenter = new G4Point3D;
#ifndef G4V10
  bbox    = new G4BoundingBox3D;
#endif
  
  G4int ii;
  for( ii = 0; ii < nPoints; ++ii ) {
    fPoints[ii] = pVec[ii];
    *cCenter += fPoints[ii];
  }
  *cCenter /= nPoints;

#ifndef G4V10
  bbox->Init(*cCenter);
  for( ii = 0; ii < nPoints; ++ii ) {
    bbox->Extend(fPoints[ii]);
  }
#endif  

  nPlanes = nFaces;
  iPlanes = new      int*[nPlanes];
  fPlanes = new G4Plane3D[nPlanes];

  G4int nn = 0;
  G4int jj, mm;
  for ( ii = 0; ii < nPlanes; ++ii ) {
    mm = theFaces[nn++];
    iPlanes[ii]    = new int[mm+1];
    iPlanes[ii][0] = mm;
    for ( jj = 0; jj < mm; jj++ )
      iPlanes[ii][jj+1] =  theFaces[nn++];
  } 

  for( ii = 0; ii < nPlanes; ii++ )
    MakePlane(ii);

//CheckConvexity(kCarTolerance/2.);
// given that we read the points from a text file, this is a more realistic precision
  CheckConvexity(0.05*micrometer);

  G4Point3D cc(0., 0., 0.);
  for( ii = 0; ii < nPoints; ii++ ) {   // recalculate the points
    MakePoint(ii);
    cc += fPoints[ii];
#ifndef G4V10
    bbox->Extend(fPoints[ii]);
#endif
  }
  cc /= nPoints;
  *cCenter = cc;

#ifndef G4V10
  bbox->Init(cc);
  for( ii = 0; ii < nPoints; ii++ ) {   // adjust bbox
    bbox->Extend(fPoints[ii]);
  }
 #endif 

  CheckConvexity(kCarTolerance/2.);   // now we can check with the default precision

  PreparePolyhedron();

}





//////////////////////////////////////////////////////////////////////////////
//
// Calculate the coef's of the plane with normal pointing to the outside
// i.e. the ThreeVectors are in anti-clockwise order when viewed from outside

void CConvexPolyhedron::MakePlane( G4int ii)
{
  G4int  np =  iPlanes[ii][0];      // how many points
  G4int *ip = &iPlanes[ii][1];      // their index into fPoints
  G4Point3D *pjj, *pkk;

// The average of the cross product of the adjacent edges of each point
// see e.g. D.F.Rogers, Proocedural Elements for Computer Graphics, p. 209
  G4Point3D fc(0., 0., 0.);       // the center of the face
  G4Point3D nv(0., 0., 0.);       // its normal vector
  for ( G4int nn = 0; nn < np; nn++ ) {
    pjj  = fPoints + ip[nn];
    pkk  = fPoints + ip[(nn+1)%np];
    nv  += G4Point3D( (pjj->y() - pkk->y()) * (pjj->z() + pkk->z()) ,
                      (pjj->z() - pkk->z()) * (pjj->x() + pkk->x()) ,
                      (pjj->x() - pkk->x()) * (pjj->y() + pkk->y()) );
    fc  += *pjj;
  }
  nv /= nv.distance();  // normalization of the normal
  fc /= np;             // the center of the face as average of its points
  G4Plane3D pl = G4Plane3D(G4Normal3D(nv), G4Point3D(fc));
  G4double  dd = pl.distance(*cCenter);
// the distance of the center of the solid must be negative to agree with the convention
// of GEANT4 which wants the normal vector to point to the outside
  if (dd > 0.) {        // adjust the equation and the order of the points
    pl = G4Plane3D(G4Normal3D(-nv), G4Point3D(fc)); 
    for ( G4int jj = 0; jj < np/2; jj++ )
      swap(ip[jj], ip[np-jj-1]);
  }
  fPlanes[ii] = pl;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate the coordinates of the point as intersection of three planes
//
//////////////////////////////////////////////////////////////////////////////
G4bool CConvexPolyhedron::MakePoint( G4int pp)
{
  G4int pl[3];
  G4int npl = 0;
  for ( G4int ii = 0; ii < nPlanes; ii++ ) {
    G4int  np =  iPlanes[ii][0];      // how many points
    G4int *ip = &iPlanes[ii][1];      // their index into fPoints
    for ( G4int jj = 0; jj < np; jj++ ) {
      if (ip[jj] == pp) {
        pl[npl++] = ii;
        break;
      }
    }
    if (npl == 3) {
       G4Point3D pnew = SolvePointEq(fPlanes[pl[0]], fPlanes[pl[1]],fPlanes[pl[2]]);
       //double dist = fPoints[pp].distance(pnew);
       fPoints[pp] = pnew;
       return true;
    }
  }
  return false; 
}

//////////////////////////////////////////////////////////////////////////////
//
// Determine the intersection of three planes
//
////////////////////////////////////////////////////////////////////////////////
G4Point3D CConvexPolyhedron::SolvePointEq(G4Plane3D& pl0, G4Plane3D& pl1, G4Plane3D& pl2)
{
  double vm[3][3], mv[3][3];
  double vd[3], pd[3];
  double det;
  G4int    ii, i1, i2, jj, j1, j2;

  vm[0][0] = pl0.a(); vm[0][1] = pl0.b(); vm[0][2] = pl0.c(); 
  vm[1][0] = pl1.a(); vm[1][1] = pl1.b(); vm[1][2] = pl1.c(); 
  vm[2][0] = pl2.a(); vm[2][1] = pl2.b(); vm[2][2] = pl2.c(); 

  vd[0] = -pl0.d(); vd[1] = -pl1.d(); vd[2] = -pl2.d(); 

  det = 0;
  ii = 0;
  i1 = 1;
  i2 = 2;
  for( jj = 0; jj < 3; jj++ ) {
    j1 = (jj+1)%3;
    j2 = (jj+2)%3;
    det += vm[ii][jj]*(vm[i1][j1]*vm[i2][j2] - vm[i1][j2]*vm[i2][j1]);
  }

  if(det == 0.) {
    G4cout << "CConvexPolyhedron::SolvePointEq : Determinant of matrix is 0." << G4endl;
    G4cerr << "CConvexPolyhedron::SolvePointEq : Determinant of matrix is 0." << G4endl;
#if defined G4V495 || G4V496 || G4V10
    G4Exception("CConvexPolyhedron::SolvePointEq","Error5",FatalException,"Determinant of matrix is 0.");
#else
    G4Exception("CConvexPolyhedron::SolvePointEq : Determinant of matrix is 0.");
#endif
  }

  det = 1. / det;
  for( ii = 0; ii < 3; ii++ ) {
    i1 = (ii+1)%3;
    i2 = (ii+2)%3;
    for( jj = 0; jj < 3; jj++ ) {
      j1 = (jj+1)%3;
      j2 = (jj+2)%3;
      mv[jj][ii] = det * (vm[i1][j1]*vm[i2][j2] - vm[i1][j2]*vm[i2][j1]);
    }
  }

#if 0
  for( ii = 0; ii < 3; ii++ ) {               // per verificare la matrice inversa
    for( jj = 0; jj < 3; jj++ ) {
      det = 0;
      for( i1 = 0; i1 < 3; i1++ ) {
        det += vm[ii][i1] * mv[i1][jj];
      }
      det = det;
    }
  }
#endif

  for( ii = 0; ii < 3; ii++ ) {
    double val = 0;
    for( jj = 0; jj < 3; jj++ ) {
      val += mv[ii][jj] * vd[jj];
    }
    pd[ii] = val;
  }
  return G4Point3D(pd[0], pd[1], pd[2]);
}


////////////////////////////////////////////////////////////////////////
//
// If the solid is convex all its defining points must not be kOutside
//
///////////////////////////////////////////////////////////////////////
void CConvexPolyhedron::CheckConvexity(const G4double dmax)
{
  for ( G4int ii =  0; ii < nPoints; ii++ ) {
    for ( G4int jj = ii; jj < nPoints; jj++ ) {
      G4ThreeVector p = G4ThreeVector((fPoints[ii] + fPoints[jj])/2.);
      for ( G4int i = 0; i < nPlanes; i++ ) {
	G4double Dist = fPlanes[i].distance(p);
	if (Dist > dmax) {
          G4cout << "WARNING - " << GetName() << " point ["
        	 << std::setw( 2) << ii << " + "
        	 << std::setw( 2) << jj << "]/2 is "
        	 << std::setiosflags(std::ios::fixed) << std::setprecision(6)
        	 << std::setw(10) << Dist/mm << " mm outside plane "
        	 << std::setw( 2) << i << G4endl;
	}
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////
//
// Shifts the plane along its normal (dist > 0 means outwards)
//
//////////////////////////////////////////////////////////////////////////
G4bool CConvexPolyhedron::MovePlane(const G4int nn, const G4double distance)
{
  if(nn < 0 || nn > nPlanes)
    return false;

  G4Plane3D pl = fPlanes[nn];
  fPlanes[nn]  = G4Plane3D(pl.a(), pl.b(), pl.c(), pl.d()-distance);
  
// needs to recalculate figure

  G4int  np =  iPlanes[nn][0];      // how many points
  G4int *ip = &iPlanes[nn][1];      // their index into fPoints
  G4Point3D fc(0., 0., 0.);       // the center of the face
  G4int ii;
  for ( ii = 0; ii < np; ii++ ) {
    MakePoint(ip[ii]);
    fc  += fPoints[ip[ii]];
  }
  fc /= np;
  //double fcDistance = fPlanes[nn].distance(fc);   // should be zero

  G4Point3D cc(0., 0., 0.);
  for( ii = 0; ii < nPoints; ii++ ) {   // recalculate the points and adjust bbox
    MakePoint(ii);
    cc += fPoints[ii];
#ifndef G4V10
    bbox->Extend(fPoints[ii]);
#endif
  }
  cc /= nPoints;
  *cCenter = cc;

  CheckConvexity(kCarTolerance/2.);   // now we can check with the default precision

  PreparePolyhedron();

  return true;
  
}

////////////////////////////////////////////////////////////////////////
//
// Prepares points and facets as used by HepPolyhedron::createPolyhedron()
//
///////////////////////////////////////////////////////////////////////////////
void CConvexPolyhedron::PreparePolyhedron()
{
// These comments taken from
// HepPolyhedron::createPolyhedron( G4int Nnodes, G4int Nfaces,
//                               const double xyz[][3],
//                               const G4int    faces[][4])
//
//***********************************************************************
//*                                                                     *
//* Name: createPolyhedron                            Date:    05.11.02 *
//* Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
//*                                                                     *
//* Function: Creates user defined polyhedron                           *
//*                                                                     *
//* Input: Nnodes  - number of nodes                                    *
//*        Nfaces  - number of faces                                    *
//*        nodes[][3] - node coordinates                                *
//*        faces[][4] - faces                                           *
//*                                                                     *
//* The faces of the polyhedron should be either triangles or planar
//* quadrilateral. Nodes of a face are defined by indexes pointing to
//* the elements in the xyz array. Numeration of the elements in the
//* array starts from 1 (like in fortran). The indexes can be positive
//* or negative. Negative sign means that the corresponding edge is
//* invisible. The normal of the face should be directed to exterior
//* of the polyhedron. 
//***********************************************************************

  G4int ii;

  if(!nFacets) {
    for (ii = 0; ii < nPlanes; ii++ ) {
      G4int *ip = iPlanes[ii];
      G4int nn = *ip;
      nFacets += 1 + (nn-3)/2;
    }
    xyz = new P3D[nPoints];
    fff = new FAC[nFacets];
  }

  for (ii = 0; ii < nPoints; ii++ ) {
    xyz[ii][0] = fPoints[ii].x();
    xyz[ii][1] = fPoints[ii].y();
    xyz[ii][2] = fPoints[ii].z();
  }

  G4int iff = 0;
  for (ii = 0; ii < nPlanes; ii++ ) {
    G4int *ip = iPlanes[ii];
    G4int nn = *ip++;
    fff[iff][0] = ip[0]+1;
    fff[iff][1] = ip[1]+1;
    fff[iff][2] = ip[2]+1;
    fff[iff][3] = (nn > 3) ? (ip[3]+1) : (0);
#if 0
    if(nn > 4) {                // to make the extra edges invisible
      fff[iff][0] *= -1;        // it works but produces a lot of annoying messages about
      fff[iff][3] *= -1;        // different edge visibility in HepPolyhedron::SetReference()
    }
#endif
    iff++;
    for ( G4int jj = 4; jj < nn; jj += 2) {
      fff[iff][0] = ip[0]+1;
      fff[iff][1] = ip[jj-1]+1;
      fff[iff][2] = ip[jj  ]+1;
      fff[iff][3] = (nn > jj+1) ? (ip[jj+1]+1) : (0);
#if 0
      fff[iff][0] *= -1;        // to make the extra edges invisible
      fff[iff][1] *= -1;
      if(nn > jj+2) {
        fff[iff][3] *= -1;
      }
#endif
      iff++;
    }
  }
}       

////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
//
//////////////////////////////////////////////////////////////////////////
G4bool CConvexPolyhedron::CalculateExtent( const EAxis pAxis,
                                           const G4VoxelLimits& pVoxelLimit,
                                           const G4AffineTransform& pTransform,
                                                 G4double& pMin, G4double& pMax) const 
{
  G4double xMin, xMax, yMin, yMax, zMin, zMax;
  G4bool flag;

  G4bool existsAfterClip = false ;
  G4ThreeVectorList*       vertices;
  pMin                   = +kInfinity;
  pMax                   = -kInfinity;
      
  vertices = CreateRotatedVertices(pTransform);   // Operator 'new' is called for vertices
      
  xMin = +kInfinity; yMin = +kInfinity; zMin = +kInfinity;
  xMax = -kInfinity; yMax = -kInfinity; zMax = -kInfinity;
      
  for( G4int nv = 0 ; nv < nPoints ; nv++ ) { 
    if( (*vertices)[nv].x() > xMax ) xMax = (*vertices)[nv].x();
    if( (*vertices)[nv].y() > yMax ) yMax = (*vertices)[nv].y();
    if( (*vertices)[nv].z() > zMax ) zMax = (*vertices)[nv].z();
      
    if( (*vertices)[nv].x() < xMin ) xMin = (*vertices)[nv].x();
    if( (*vertices)[nv].y() < yMin ) yMin = (*vertices)[nv].y();
    if( (*vertices)[nv].z() < zMin ) zMin = (*vertices)[nv].z();
  }

  if ( pVoxelLimit.IsXLimited() ) {
    if ( (xMin > pVoxelLimit.GetMaxXExtent() + kCarTolerance)
      || (xMax < pVoxelLimit.GetMinXExtent() - kCarTolerance) ) {
      delete vertices ;    //  'new' in the function called
      return false ;
    } 
    else {
      if ( xMin < pVoxelLimit.GetMinXExtent() ) {
        xMin = pVoxelLimit.GetMinXExtent() ;
      }
      if ( xMax > pVoxelLimit.GetMaxXExtent() ) {
        xMax = pVoxelLimit.GetMaxXExtent() ;
      }
    }
  }

  if ( pVoxelLimit.IsYLimited() ) {
    if ( (yMin > pVoxelLimit.GetMaxYExtent() + kCarTolerance)
      || (yMax < pVoxelLimit.GetMinYExtent() - kCarTolerance) ) {
      delete vertices ;    //  'new' in the function called
      return false;
    }
    else {
      if ( yMin < pVoxelLimit.GetMinYExtent() ) {
        yMin = pVoxelLimit.GetMinYExtent() ;
      }
      if ( yMax > pVoxelLimit.GetMaxYExtent() ) {
        yMax = pVoxelLimit.GetMaxYExtent() ;
      }
    }
  }

  if ( pVoxelLimit.IsZLimited() ) {
    if ( (zMin > pVoxelLimit.GetMaxZExtent() + kCarTolerance)
      || (zMax < pVoxelLimit.GetMinZExtent() - kCarTolerance) ) {
      delete vertices ;    //  'new' in the function called
      return false;
    }
    else {
      if ( zMin < pVoxelLimit.GetMinZExtent() ) {
        zMin = pVoxelLimit.GetMinZExtent() ;
      }
      if ( zMax > pVoxelLimit.GetMaxZExtent() ) {
        zMax = pVoxelLimit.GetMaxZExtent() ;
      }
    }
  } 

  switch (pAxis)
  {
    case kXAxis:
    pMin=xMin;
    pMax=xMax;
    break;

    case kYAxis:
      pMin=yMin;
      pMax=yMax;
      break;

    case kZAxis:
      pMin=zMin;
      pMax=zMax;
      break;

    default:
      break;
  }

  if ( (pMin != kInfinity) || (pMax != -kInfinity) ) {
    existsAfterClip=true;
    // Add tolerance to avoid precision troubles
    pMin -= kCarTolerance ;
    pMax += kCarTolerance ;      
  }

  delete vertices ;          //  'new' was called in CreateRotatedVertices()
  flag = existsAfterClip ;

  return flag;
}

////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance
//
////////////////////////////////////////////////////////////////////////
EInside CConvexPolyhedron::Inside( const G4ThreeVector& p ) const 
{
  EInside in = kInside;
  G4double Dist;
  G4int i;

#ifndef G4V10
  i = bbox->Inside(p);
  if (!i)
    return in = kOutside;
#endif

  for ( i = 0; i < nPlanes; i++ ) {
    Dist = fPlanes[i].distance(p);
    if (Dist > kCarTolerance/2) {
      return in=kOutside;
    }
    else if (Dist>-kCarTolerance/2) {
       in=kSurface;
    } 
  }
  return in;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If 2+ sides equidistant, first side's normal returned (arbitrarily)
//
///////////////////////////////////////////////////////////////////////////////
G4ThreeVector CConvexPolyhedron::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4double safe = kInfinity;
  G4double Dist;
  G4int    i, imin = 0;

  for ( i = 0; i < nPlanes; i++ ) {
    Dist = fabs(fPlanes[i].distance(p));
    if (Dist < safe) {
      safe=Dist;
      imin=i;
    }
  }
  return G4ThreeVector(fPlanes[imin].normal());
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside - return kInfinity if no intersection
//
// ALGORITHM:
// For each component, calculate pair of minimum and maximum intersection
// values for which the particle is in the extent of the shape
// - The smallest (MAX minimum) allowed distance of the pairs is intersect
//
//////////////////////////////////////////////////////////////////////////////////
G4double CConvexPolyhedron::DistanceToIn( const G4ThreeVector& p,  const G4ThreeVector& v ) const 
{
  G4double snxt;    // snxt = default return value
  G4double smax, smin;
  G4double pdist, Comp, vdist;

  smin = 0;
  smax = kInfinity;

  for ( G4int i = 0; i < nPlanes; i++ ) {
    pdist = fPlanes[i].distance(p);
    Comp  = fPlanes[i].a()*v.x() + fPlanes[i].b()*v.y() + fPlanes[i].c()*v.z();
    if (pdist >= -0.5*kCarTolerance) {  // Outside the plane -> this is an extent entry distance
      if (Comp >= 0) {
        return snxt = kInfinity;
      }
      else {
        vdist = -pdist / Comp;
        if (vdist > smin) {
          if (vdist < smax)
            smin = vdist;
          else
            return snxt = kInfinity;
        }
      }
    }
    else {                              // Inside the plane -> coud be an extent exit distance (smax)
      if (Comp > 0) {                   // Will leave extent
        vdist = -pdist / Comp;
        if (vdist < smax) {
          if (vdist > smin)
            smax = vdist;
          else
            return snxt = kInfinity;
        }  
      }
    }
  }

  if (smin >= 0 ) {
    snxt = smin;
  }
  else {
    snxt = 0;
  }

#if 0
  G4ThreeVector mid = p + 0.5*(smax+smin)*v;
  if(Inside(mid) == kOutside) {
    // cannot happen if solid is convex
    G4cout << " CConvexPolyhedron::DistanceToIn(p,v) middle point outside" << G4endl;
    snxt=kInfinity;
  }
#endif

  return snxt;

}

///////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from outside
// This is the best fast estimation of the shortest distance to trap
// - Returns 0 is ThreeVector inside
//
//////////////////////////////////////////////////////////////////////////////
G4double CConvexPolyhedron::DistanceToIn( const G4ThreeVector& p ) const 
{
  G4double safe,Dist;
  G4int i;

  safe = kInfinity;
  for ( i = 0; i < nPlanes; i++ ) {
    Dist = fPlanes[i].distance(p);
    if (Dist > -kCarTolerance/2 && Dist < safe)
      safe = Dist;
  }
  if (safe < 0) safe=0;
  return safe; 
}

/////////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from inside
// Calculate distance to x/y/z planes - smallest is exiting distance
//
//////////////////////////////////////////////////////////////////////////////////////
G4double CConvexPolyhedron::DistanceToOut( const G4ThreeVector& p,  const G4ThreeVector& v,
                                           const G4bool calcNorm,
                                                 G4bool *validNorm, G4ThreeVector *n) const 
{
  G4double snxt = kInfinity;
  G4double pdist,Comp,vdist;
  G4int    i, side = -1;

  for ( i = 0; i < nPlanes; i++ ) {
    pdist = fPlanes[i].distance(p);
    Comp  = fPlanes[i].a()*v.x()+fPlanes[i].b()*v.y()+fPlanes[i].c()*v.z();
    if (pdist > 0) {                        // Outside the plane
      if (Comp > 0) {                       // Leaving immediately
        if (calcNorm) {
          *validNorm = true;
          *n = G4ThreeVector(fPlanes[i].a(), fPlanes[i].b(), fPlanes[i].c());
        }
        return snxt = 0;
      }
    }
    else if (pdist < -0.5*kCarTolerance) {  // Inside the plane
      if (Comp > 0) {                       // Will leave extent
        vdist = -pdist / Comp;
        if (vdist < snxt) {
          snxt = vdist;
          side = i;
        }
      }
    }
    else {                                  // On surface
      if (Comp > 0) {
        if (calcNorm) {
          *validNorm = true;
          *n = G4ThreeVector(fPlanes[i].a(), fPlanes[i].b(), fPlanes[i].c());
        }
        return snxt = 0;
      }
    }
  }

  if (calcNorm) {
    *validNorm = true;
    *n = G4ThreeVector(fPlanes[side].a(), fPlanes[side].b(), fPlanes[side].c());
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is ThreeVector outside
//
////////////////////////////////////////////////////////////////////////////////
G4double CConvexPolyhedron::DistanceToOut( const G4ThreeVector& p ) const 
{
  G4double safe,Dist;
  G4int i;

#ifdef G4CSGDEBUG
  if( Inside(p) == kOutside )
  {
     G4cout.precision(16) ;
     G4cout << G4endl ;
     DumpInfo();
     G4cout << "Position:"  << G4endl << G4endl ;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
     G4cout << "CConvexPolyhedron::DistanceToOut(p) - point p is outside ?!" << G4endl ;
     G4cerr << "CConvexPolyhedron::DistanceToOut(p) - point p is outside ?!" << G4endl ;
  }
#endif

  safe = kInfinity;
  for( i = 0; i < nPlanes; i++ ) {
    Dist = -fPlanes[i].distance(p);
    if(Dist > -kCarTolerance/2 && Dist < safe)
      safe = Dist;
  }
  if (safe < 0) safe=0;
  return safe;
}

//////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Note: Caller has deletion resposibility
//
///////////////////////////////////////////////////////////////////////////////
G4ThreeVectorList*
CConvexPolyhedron::CreateRotatedVertices( const G4AffineTransform& pTransform ) const 
{
  G4ThreeVectorList *vertices;
  vertices = new G4ThreeVectorList();
  vertices->reserve(nPoints);
  if (vertices) {
    for ( G4int i = 0; i < nPoints; i++ ) {
      vertices->push_back(pTransform.TransformPoint(fPoints[i]));
    }
  }
  else {
//  DumpInfo();
#if defined G4V495 || G4V496 || G4V10
    G4Exception("CConvexPolyhedron::CreateRotatedVertices()","error6", FatalException,"Out of memory !");
#else
    G4Exception("CConvexPolyhedron::CreateRotatedVertices() - Out of memory !");
#endif
  }
  return vertices;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType
//
//////////////////////////////////////////////////////////////////////////////
G4GeometryType CConvexPolyhedron::GetEntityType() const
{
  return G4String("CConvexPolyhedron");
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream
//
//////////////////////////////////////////////////////////////////////////////
std::ostream& CConvexPolyhedron::StreamInfo( std::ostream& os ) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: CConvexPolyhedron\n"
     << nPoints << " Points: \n"
     << nPlanes << " Planes: \n"
     << "-----------------------------------------------------------\n";

  return os;
}

///////////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation
//
//////////////////////////////////////////////////////////////////////////////
void CConvexPolyhedron::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
#ifdef G4V47
  scene.AddSolid (*this);
#else
  //  scene.AddThis (*this);
  scene.AddSolid (*this);
#endif
}

G4Polyhedron* CConvexPolyhedron::CreatePolyhedron () const
{
  G4Polyhedron *cpp = new G4Polyhedron();
  cpp->createPolyhedron(nPoints, nFacets, xyz, fff);
  return cpp;
}

#ifndef G4V10
G4NURBS* CConvexPolyhedron::CreateNURBS () const
{
  return 0;
}
#endif
