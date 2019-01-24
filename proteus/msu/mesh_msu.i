  /* mesh_msu.i */

deprecated

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