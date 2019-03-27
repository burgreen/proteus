// msu_mesh_util.h

#ifndef msu_mesh_util_h
#define msu_mesh_util_h

#include <string>

#include "msu/msu_core.h"

#include "mesh.h"

int msu_mesh_util_find_index_elementBoundaryMaterialTypes( int&             idx,
                                                           Mesh&            mesh, 
                                                           const msu_vecI&  vecN,
                                                           std::string      cmd="undefined",
                                                           double           tol=1.e-6 );

int msu_mesh_util_find_index_elementBoundaryMaterialTypes( int&             idx,
                                                           Mesh&            mesh, 
                                                           const msu_vecD&  centroid,
                                                           std::string      cmd="undefined",
                                                           double           tol=1.e-5 );
int msu_mesh_util_area_normal( Mesh&            mesh, 
                               const msu_vecI&  vecN,
                               double&          area,
                               double*          normal );

int msu_mesh_util_area_normal_2d( Mesh&            mesh, 
                                  const msu_vecI&  vecN,
                                  double&          area,
                                  double*          normal );

int msu_mesh_util_area_normal( msu_vecD&        x, 
                               msu_vecD&        y, 
                               msu_vecD&        z, 
                               const msu_vecI&  vecN,
                               double&          area,
                               double*          normal );

int msu_mesh_util_area_normal_2d( msu_vecD&        x, 
                                  msu_vecD&        y, 
                                  msu_vecD&        z, 
                                  const msu_vecI&  vecN,
                                  double&          area,
                                  double*          normal );

#endif
