from libcpp cimport bool

cimport mesh as cppm

cdef extern from "msu/mesh_msu.h":

  cdef int  mesh_msu_readBC( cppm.Mesh& mesh, const char* filebase, int indexBase );

  cdef int  mesh_msu_readBC_tetgen( cppm.Mesh&  mesh, 
                                    const char* filebase, 
                                    int         indexBase );
  
  cdef int  mesh_msu_readBC_tetgen_parallel( cppm.Mesh&  mesh, 
                                             const char* filebase, 
                                             int         indexBase );
  
  cdef int  mesh_msu_readBC_triangle( cppm.Mesh&  mesh, 
                                      const char* filebase, 
                                      int         num_interiorEdgeTags, 
                                      int*        interiorEdgeTags );

  cdef int  mesh_msu_readBC_hex( cppm.Mesh&  mesh, 
                                 const char* filebase, 
                                 int         indexBase );

  cdef int  mesh_msu_read_hex( cppm.Mesh&  mesh, 
                               const char* filebase, 
                               int         indexBase );
  