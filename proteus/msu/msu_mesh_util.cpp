// msu_mesh_util.cpp

#include <ANN/ANN.h>

#include "msu/msu_core.h"

#include "mesh.h"

using namespace std;

static ANNkd_tree* s_tree;

#include "msu_mesh_util.ipp"

static int s_interior_bc = 0; // hardwire

  int 
//---------------------------------------------------------------
  msu_mesh_util_find_index_elementBoundaryMaterialTypes
//---------------------------------------------------------------
(
  int&             idx,
  Mesh&            mesh, 
  const msu_vecD&  centroid,
  std::string      cmd,
  double           tol
)
{
  int nbF = mesh.nExteriorElementBoundaries_global;
  int niF = mesh.nInteriorElementBoundaries_global;

  //- - - - - - - - - - - - - - - - - - - - - - - 
  if( cmd == "init" )
  {
    int nN  = mesh.nNodes_elementBoundary;
    int nF  = nbF + niF;
    
    FOR( f, nbF )
    {
      int idx = mesh.exteriorElementBoundariesArray[f];
      mesh.elementBoundaryMaterialTypes[idx] = s_interior_bc;
    }

    FOR( f, niF )
    {
      int idx = mesh.interiorElementBoundariesArray[f];
      mesh.elementBoundaryMaterialTypes[idx] = s_interior_bc;
    }

std::cout << "msu_mesh_util: nN  = " << nN << std::endl;
std::cout << "msu_mesh_util: nbF = " << nbF << std::endl;
std::cout << "msu_mesh_util: niF = " << niF << std::endl;
std::cout << "msu_mesh_util: nF  = " << nF << std::endl;

    ANNpointArray centroids = annAllocPts( nF, 3 );

    double centroid[3]; 

    FOR( f, nbF )
    {
      int eF = mesh.exteriorElementBoundariesArray[f];
      FOR(i,3) centroid[i] = 0;
      FOR(j,nN)
      {
        int n = mesh.elementBoundaryNodesArray[eF*nN+j];
        FOR(i,3) centroid[i] += mesh.nodeArray[n*3+i] / nN;
      }
      centroids[f][0] = centroid[0];
      centroids[f][1] = centroid[1];
      centroids[f][2] = centroid[2];
    }

    FOR( f, niF )
    {
      int iF = mesh.interiorElementBoundariesArray[f];
      FOR(i,3) centroid[i] = 0;
      FOR(j,nN)
      {
        int n = mesh.elementBoundaryNodesArray[iF*nN+j];
        FOR(i,3) centroid[i] += mesh.nodeArray[n*3+i] / nN;
      }
      centroids[nbF+f][0] = centroid[0];
      centroids[nbF+f][1] = centroid[1];
      centroids[nbF+f][2] = centroid[2];
    }
  
    s_tree = new ANNkd_tree( centroids, nF, 3 );

    return msu_success;
  }

  //- - - - - - - - - - - - - - - - - - - - - - - 
  if( cmd == "done" )
  {
    if( s_tree ) delete s_tree; 
    annClose(); 
    return msu_success;
  }

  //- - - - - - - - - - - - - - - - - - - - - - - 
  // cmd == "find"  locate the idx of the centroid

  idx = -1;

  int find_at_least = 1; 
  int find_at_most = 1; 
  ANNidxArray nnIdx = new ANNidx[find_at_most];
  ANNdistArray nnDists = new ANNdist[find_at_most];
  ANNpoint pt = annAllocPt(3);

  FOR(i,3) pt[i] = centroid[i];

  double r   = tol; // srch radius 
  double r2  = r * r;
  double eps = 0;

  int found = s_tree->annkFRSearch( pt, r2, find_at_most, nnIdx, nnDists, eps );

/* parallel version
*/
  if( found < find_at_least )
  {
    //std::cerr << "found < find_at_least" << std::endl;
    //return msu_failure;
  }

  if( found > find_at_most )
  {
    std::cerr << "found > find_at_most" << std::endl;
    return msu_failure;
  }

  if( found == 1 )
  {
    int f = nnIdx[0];

    if( f < nbF ) idx = mesh.exteriorElementBoundariesArray[f];
    else          idx = mesh.interiorElementBoundariesArray[f-nbF];
  }

  delete [] nnIdx;
  delete [] nnDists;
  
  return msu_success;
}

  int 
//---------------------------------------------------------------
  msu_mesh_util_find_index_elementBoundaryMaterialTypes
//---------------------------------------------------------------
(
  int&             idx,
  Mesh&            mesh, 
  const msu_vecI&  vecN,
  std::string      cmd,
  double           tol
)
{
  int nbF = mesh.nExteriorElementBoundaries_global;
  int niF = mesh.nInteriorElementBoundaries_global;

  //----------------------------------------------
  if( cmd == "init" )
  {
    int nN  = mesh.nNodes_elementBoundary;
    int nF  = nbF + niF;
    
std::cout << "msu_mesh_util: nN  = " << nN << std::endl;
std::cout << "msu_mesh_util: nbF = " << nbF << std::endl;
std::cout << "msu_mesh_util: niF = " << niF << std::endl;
std::cout << "msu_mesh_util: nF  = " << nF << std::endl;

    ANNpointArray centroids = annAllocPts( nF, 3 );

    double centroid[3]; 

    FOR( f, nbF )
    {
      int eF = mesh.exteriorElementBoundariesArray[f];
      FOR(i,3) centroid[i] = 0;
      FOR(j,nN)
      {
        int n = mesh.elementBoundaryNodesArray[eF*nN+j];
        FOR(i,3) centroid[i] += mesh.nodeArray[n*3+i] / nN;
      }
      centroids[f][0] = centroid[0];
      centroids[f][1] = centroid[1];
      centroids[f][2] = centroid[2];
    }

    FOR( f, niF )
    {
      int iF = mesh.interiorElementBoundariesArray[f];
      FOR(i,3) centroid[i] = 0;
      FOR(j,nN)
      {
        int n = mesh.elementBoundaryNodesArray[iF*nN+j];
        FOR(i,3) centroid[i] += mesh.nodeArray[n*3+i] / nN;
      }
      centroids[nbF+f][0] = centroid[0];
      centroids[nbF+f][1] = centroid[1];
      centroids[nbF+f][2] = centroid[2];
    }
  
    s_tree = new ANNkd_tree( centroids, nF, 3 );

    FOR( f, nF )  mesh.elementBoundaryMaterialTypes[f] = s_interior_bc;

    return msu_success;
  }

  //----------------------------------------------
  if( cmd == "done" )
  {
    if( s_tree ) delete s_tree; 
    annClose(); 
    return msu_success;
  }

  //----------------------------------------------
  // cmd == "find"  locate the idx of the centroid

  idx = -1;

  msu_vecD centroid(3,0); 
  {
    FOR( j, vecN.size() ) FOR(i,3) centroid[i] += mesh.nodeArray[vecN[j]*3+i] / vecN.size();
  }

  int find_at_least = 1; 
  int find_at_most = 1; 
  ANNidxArray nnIdx = new ANNidx[find_at_most];
  ANNdistArray nnDists = new ANNdist[find_at_most];
  ANNpoint pt = annAllocPt(3);

  FOR(i,3) pt[i] = centroid[i];

  double r   = tol; // srch radius 
  double r2  = r * r;
  double eps = 0;

  int found = s_tree->annkFRSearch( pt, r2, find_at_most, nnIdx, nnDists, eps );

/* parallel version
*/
  if( found < find_at_least )
  {
    //std::cerr << "found < find_at_least" << std::endl;
    //return msu_failure;
  }

  if( found > find_at_most )
  {
    std::cerr << "found > find_at_most" << std::endl;
    return msu_failure;
  }

  if( found == 1 )
  {
    int f = nnIdx[0];

    if( f < nbF ) idx = mesh.exteriorElementBoundariesArray[f];
    else          idx = mesh.interiorElementBoundariesArray[f-nbF];
  }

  delete [] nnIdx;
  delete [] nnDists;
  
  return msu_success;
}

  int 
//---------------------------------------------------------------
  msu_mesh_util_area_normal
//---------------------------------------------------------------
(
  Mesh&            mesh, 
  const msu_vecI&  vecN,
  double&          area,
  double*          normal
)
{
  s_face_area_normal( vecN, mesh.nodeArray, area, normal );

  return msu_success;
}

  int 
//---------------------------------------------------------------
  msu_mesh_util_area_normal
//---------------------------------------------------------------
(
  msu_vecD&        x, 
  msu_vecD&        y, 
  msu_vecD&        z, 
  const msu_vecI&  vecN,
  double&          area,
  double*          normal
)
{
  s_face_area_normal( vecN, x,y,z, area, normal );

  return msu_success;
}

  int 
//---------------------------------------------------------------
  msu_mesh_util_area_normal_2d
//---------------------------------------------------------------
(
  Mesh&            mesh, 
  const msu_vecI&  vecN,
  double&          area,
  double*          normal
)
{
  if( vecN.size() != 2 ) return msu_failure;

  double tol = 1.e-4; // hardwire

  int n0 = vecN[0];
  int n1 = vecN[1];

  double x0 = mesh.nodeArray[n0*3+0];
  double y0 = mesh.nodeArray[n0*3+1];
  double z0 = mesh.nodeArray[n0*3+2];

  double x1 = mesh.nodeArray[n1*3+0];
  double y1 = mesh.nodeArray[n1*3+1];
  double z1 = mesh.nodeArray[n1*3+2];

  double dx = x1 - x0;
  double dy = y1 - y0;

  area = std::sqrt( dx*dx + dy*dy );
  normal[0] = -dy/area;
  normal[1] = +dx/area;
  normal[2] = 0;

  return msu_success;
}

  int 
//---------------------------------------------------------------
  msu_mesh_util_area_normal_2d
//---------------------------------------------------------------
(
  msu_vecD&        x, 
  msu_vecD&        y, 
  msu_vecD&        z, 
  const msu_vecI&  vecN,
  double&          area,
  double*          normal
)
{
  if( vecN.size() != 2 ) return msu_failure;

  double tol = 1.e-4; // hardwire

  int n0 = vecN[0];
  int n1 = vecN[1];

  double x0 = x[n0];
  double y0 = y[n0];
  double z0 = z[n0];

  double x1 = x[n1];
  double y1 = y[n1];
  double z1 = z[n1];

  double dx = x1 - x0;
  double dy = y1 - y0;

  area = std::sqrt( dx*dx + dy*dy );
  normal[0] = -dy/area;
  normal[1] = +dx/area;
  normal[2] = 0;

  return msu_success;
}
