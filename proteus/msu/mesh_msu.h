  /* mesh_msu.h */

#ifndef msu_mesh_msu_h
#define msu_mesh_msu_h

  struct Mesh;

  int  mesh_msu_readBC( Mesh& mesh, const char* filebase, int indexBase );

  int  mesh_msu_readBC_tetgen( Mesh&       mesh, 
                               const char* filebase, 
                               int         indexBase );
  
  int  mesh_msu_readBC_tetgen_parallel( Mesh&       mesh, 
                                        const char* filebase, 
                                        int         indexBase );
  
  int  mesh_msu_readBC_triangle( Mesh&       mesh, 
                                 const char* filebase, 
                                 int         num_interiorEdgeTags, 
                                 int*        interiorEdgeTags );

  int  mesh_msu_readBC_hex( Mesh&       mesh, 
                            const char* filebase, 
                            int         indexBase );

  int  mesh_msu_read_hex( Mesh& mesh, const char* filebase, int indexBase );

  extern int g_mesh_msu;

#endif
