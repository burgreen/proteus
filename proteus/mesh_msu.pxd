from libcpp cimport bool

cimport mesh as cppm

cdef extern from "msu/mesh_msu.h":

  cdef cppm.Mesh Mesh

  cdef int  mesh_msu_readBC( Mesh& mesh, const char* filebase, int indexBase );

  cdef int  mesh_msu_readBC_tetgen( Mesh&       mesh, 
                                    const char* filebase, 
                                    int         indexBase );
  
  cdef int  mesh_msu_readBC_tetgen_parallel( Mesh&       mesh, 
                                             const char* filebase, 
                                             int         indexBase );
  
  cdef int  mesh_msu_readBC_triangle( Mesh&       mesh, 
                                      const char* filebase, 
                                      int         num_interiorEdgeTags, 
                                      int*        interiorEdgeTags );

  cdef int  mesh_msu_readBC_hex( Mesh&       mesh, 
                                 const char* filebase, 
                                 int         indexBase );

  cdef int  mesh_msu_read_hex( Mesh&       mesh, 
                               const char* filebase, 
                               int         indexBase );
  