// mesh_msu.cpp

#include <map>
#include <set>
#include <iostream>

#include "mesh.h"

#include "msu/msu_core.h"
#include "msu/msu_mesh_util.h"

static std::istream& s_eatline( std::istream& s )
{
  while (s.get() != '\n' && s.good()) {}
  return s;
}

static std::istream& s_eatcomments( std::istream& s ) 
{
  char c = s.peek() ;
  //mwf add whitespace too, or just blank lines?
  while ( ( c == '\n' || c == '!' || c == '%' || c == '#' || c == ';' || c == '$')
    && s.good() ) 
    { s >> s_eatline ; c = s.peek() ; }
  return s ;
}

static int s_interior_bc = 0; // hardwire

  int 
//---------------------------------------------------------------
  mesh_msu_readBC
//---------------------------------------------------------------
// read node and element BC information 
(
  Mesh&        mesh, 
  const char*  filebase, 
  int          indexBase
)
{
  using namespace std;

  assert(filebase);

  bool failed = true;

  std::string bcFilename = std::string(filebase)+".bc";

  std::ifstream bcFile( bcFilename.c_str() );
  if( ! bcFile.good() )
  {
    std::cerr<<"readBC cannot open file " << bcFilename << std::endl;
    failed = true;
    return failed;
  }

  std::string firstWord;
  bcFile>>firstWord;
  if( firstWord != "BC" )
  {
    std::cerr<<"firstWord is not BC" << std::endl;
    failed = true;
    return failed;
  }

  int numN, numT, numQ, bc;
  bcFile >> numN >> numT >> numQ;

  int n, index; std::vector<int> vecN;

  failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN, "init" );
  if( failed ) { std::cerr << "bad init" << std::endl; return failed; }

  // 1. boundary nodes

  FOR( i, numN )
  {
    bcFile >> n >> bc;
    mesh.nodeMaterialTypes[n-indexBase] = bc;
  }

  // 2. boundary triangles

  vecN.resize(3);

  FOR( i, numT )
  {
    bcFile >> vecN[0] >> vecN[1] >> vecN[2] >> bc;
    for( int j=0; j<3; j++ ) vecN[j] -= indexBase;
    failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN );
    if( failed ) { std::cerr << "did not find tri bnd face index" << std::endl; return failed; }
    //mesh.elementBoundaryMaterialTypes[index] = flag;
    mesh.elementBoundaryMaterialTypes[index] = bc;
  }

  // 3. boundary quads

  vecN.resize(4);

  FOR( i, numQ )
  {
    bcFile >> vecN[0] >> vecN[1] >> vecN[2] >> vecN[3] >> bc;
    for( int j=0; j<4; j++ ) vecN[j] -= indexBase;
    failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN );
    if( failed ) { std::cerr << "did not find quad bnd face index" << std::endl; return failed; }
    //mesh.elementBoundaryMaterialTypes[index] = flag;
    mesh.elementBoundaryMaterialTypes[index] = bc;
  }

  msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN, "done" );

  return msu_success;
}

  int 
//---------------------------------------------------------------
  mesh_msu_readBC_tetgen
//---------------------------------------------------------------
// read node and element BC information 
(
  Mesh&        mesh, 
  const char*  filebase, 
  int          indexBase
)
{
  using namespace std;

  assert(filebase);

  bool failed = true;

  std::string bcFilename = std::string(filebase)+".face";

  std::ifstream bcFile( bcFilename.c_str() );

  if( ! bcFile.good() )
  {
    std::cerr<<"readBC cannot open file " << bcFilename << std::endl;
    failed = true;
    return failed;
  }

  int numF, numA;
  bcFile >> s_eatcomments >> numF >> numA;

  int n, index; std::vector<int> vecN;

  failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN, "init" );

  if( failed ) { std::cerr << "bad init" << std::endl; return failed; }

  // 1. boundary triangles

  int idx, bc;
  double area, normal[3];

  vecN.resize(3);

  std::set<int> attrF;
  std::map<int,msu_vecD> nx;
  std::map<int,msu_vecD> ny;
  std::map<int,msu_vecD> nz;

  attrF.insert(s_interior_bc);

  int numF_assigned = 0;

  FOR( i, numF )
  {
    // a. set bc 

    bcFile >> idx >> vecN[0] >> vecN[1] >> vecN[2] >> bc;

    for( int j=0; j<3; j++ ) vecN[j] -= indexBase;
    failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN );
    if( failed ) { std::cerr << "did not locate bnd face centroid" << std::endl; return failed; }
    if( index > -1 )
    {
      mesh.elementBoundaryMaterialTypes[index] = bc;
      numF_assigned++;
    }

    attrF.insert(bc);

    // b. compute face normal

    msu_mesh_util_area_normal( mesh, vecN, area, normal );

    nx[bc].push_back( normal[0] );
    ny[bc].push_back( normal[1] );
    nz[bc].push_back( normal[2] );
  }

  std::cout << "serial numF_assigned = " << numF_assigned << std::endl;

  // compute mean and stddev of bc zone face normals

  int nA = attrF.size();
  mesh.num_attrF = nA;
  mesh.attrF     = new int[nA];
  mesh.nx_mean   = new double[nA];
  mesh.ny_mean   = new double[nA];
  mesh.nz_mean   = new double[nA];
  mesh.nx_stddev = new double[nA];
  mesh.ny_stddev = new double[nA];
  mesh.nz_stddev = new double[nA];

  int cnt = 0;
  std::set<int>::iterator I;
  for( I = attrF.begin(); I != attrF.end(); ++I )
  {
    int a = *I;
    //std::cout << "attrF = " << a << std::endl;

    double nx_mean=0, nx_stddev=0;
    double ny_mean=0, ny_stddev=0;
    double nz_mean=0, nz_stddev=0;

    if( a != s_interior_bc )
    {
      FOR( i, nx[a].size() ) nx_mean += nx[a][i];
      FOR( i, ny[a].size() ) ny_mean += ny[a][i];
      FOR( i, nz[a].size() ) nz_mean += nz[a][i];
  
      nx_mean /= nx[a].size();
      ny_mean /= ny[a].size();
      nz_mean /= nz[a].size();
  
      FOR( i, nx[a].size() ) nx_stddev += ((nx[a][i] - nx_mean) * (nx[a][i] - nx_mean));
      FOR( i, ny[a].size() ) ny_stddev += ((ny[a][i] - ny_mean) * (ny[a][i] - ny_mean));
      FOR( i, nz[a].size() ) nz_stddev += ((nz[a][i] - nz_mean) * (nz[a][i] - nz_mean));
  
      nx_stddev = std::sqrt( nx_stddev / (nx[a].size()-1) );
      ny_stddev = std::sqrt( ny_stddev / (ny[a].size()-1) );
      nz_stddev = std::sqrt( nz_stddev / (nz[a].size()-1) );
    }

    mesh.attrF[cnt]     = a;
    mesh.nx_mean[cnt]   = nx_mean;
    mesh.ny_mean[cnt]   = ny_mean;
    mesh.nz_mean[cnt]   = nz_mean;
    mesh.nx_stddev[cnt] = nx_stddev;
    mesh.ny_stddev[cnt] = ny_stddev;
    mesh.nz_stddev[cnt] = nz_stddev;
    cnt++;
  }

  msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN, "done" );

  return msu_success;
}

  int 
//---------------------------------------------------------------
  mesh_msu_readBC_tetgen_parallel
//---------------------------------------------------------------
// read node and element BC information 
(
  Mesh&           mesh, 
  const char*     filebase, 
  int             indexBase
)
{
  using namespace std;

  assert(filebase);

  bool failed = true;

  std::string bcFilename = std::string(filebase)+".bcmsh";

  std::ifstream ifs( bcFilename.c_str() );

  if( ! ifs.good() )
  {
    std::cerr<<"readBC cannot open file " << bcFilename << std::endl;
    failed = true;
    return failed;
  }

  // 0. initialize the local mesh partition

  int n, index; std::vector<double> centroid(3);

  failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, centroid, "init" );

  if( failed ) { std::cerr << "bad init" << std::endl; return failed; }

  // 1. find each triangle defined in .bcmsh file

  int nN;
  double area, normal[3];

  vector<int> vecN(3);

  std::set<int> attrF;
  std::map<int,msu_vecD> nx;
  std::map<int,msu_vecD> ny;
  std::map<int,msu_vecD> nz;

  attrF.insert(s_interior_bc);

  int numF_assigned = 0;

  // read bcmsh file

  std::string text;
  int version; ifs >> text >> version; // std::cout << text << version << std::endl;
  int dim; ifs >> text >> dim;         // std::cout << text << dim << std::endl;
  int numn; ifs >> text >> numn;       // std::cout << text << numn << std::endl;
  int numf; ifs >> text >> numf;          std::cout << text << numf << std::endl;
  int numN; ifs >> text >> numN;       // std::cout << text << numN << std::endl;
  msu_vecD x(numN);
  msu_vecD y(numN);
  msu_vecD z(numN);
  FOR( i, numN ) ifs >> x[i] >> y[i] >> z[i];
  int numA; ifs >> text >> numA;       // std::cout << text << numA << std::endl;
  msu_vecI vec_attrF(numA);
  FOR( i, numA ) ifs >> vec_attrF[i];
  int ordered; ifs >> text >> ordered; // std::cout << text << ordered << std::endl;
  int numF; ifs >> text >> numF;       // std::cout << text << numF << std::endl;

  FOR( i, numF )
  {
    // a. set bc 

    ifs >> nN >> vecN[0] >> vecN[1] >> vecN[2];
    FOR(j,3) vecN[j]--; // base0
    FOR(j,3) centroid[j] = 0;
    FOR(j,3) centroid[0] += x[vecN[j]]; 
    FOR(j,3) centroid[1] += y[vecN[j]]; 
    FOR(j,3) centroid[2] += z[vecN[j]]; 
    FOR(j,3) centroid[j] /= 3;

    int bc = vec_attrF[i];

    failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, centroid );
    if( failed ) { std::cerr << "did not locate bnd face centroid" << std::endl; return failed; }

    if( index > -1 ) // therefore, face is in the mesh partition
    {
      mesh.elementBoundaryMaterialTypes[index] = bc;

      attrF.insert(bc);

      // b. compute face normal

      msu_mesh_util_area_normal( x,y,z, vecN, area, normal );

      nx[bc].push_back( normal[0] );
      ny[bc].push_back( normal[1] );
      nz[bc].push_back( normal[2] );

      numF_assigned++;
    }
  }

  std::cout << "parallel numF_assigned = " << numF_assigned << std::endl;
  std::map<int,int> mapIC; // idx-count
  int nbF = mesh.nExteriorElementBoundaries_global;
  FOR( f, nbF )
  {
    int idx = mesh.exteriorElementBoundariesArray[f];
    int bc = mesh.elementBoundaryMaterialTypes[idx];
    mapIC[bc]++;
  }
  std::map<int,int>::iterator IC; 
  ITER( IC, mapIC ) std::cout << "mapIC " << IC->first <<": "<< IC->second << std::endl;

  // compute mean and stddev of bc zone face normals

  int nA = attrF.size();
  mesh.num_attrF = nA;
  mesh.attrF     = new int[nA];
  mesh.nx_mean   = new double[nA];
  mesh.ny_mean   = new double[nA];
  mesh.nz_mean   = new double[nA];
  mesh.nx_stddev = new double[nA];
  mesh.ny_stddev = new double[nA];
  mesh.nz_stddev = new double[nA];

  int cnt = 0;
  std::set<int>::iterator I;
  for( I = attrF.begin(); I != attrF.end(); ++I )
  {
    int a = *I;

    double nx_mean=0, nx_stddev=0;
    double ny_mean=0, ny_stddev=0;
    double nz_mean=0, nz_stddev=0;

    if( a != s_interior_bc )
    {
      FOR( i, nx[a].size() ) nx_mean += nx[a][i];
      FOR( i, ny[a].size() ) ny_mean += ny[a][i];
      FOR( i, nz[a].size() ) nz_mean += nz[a][i];
  
      nx_mean /= nx[a].size();
      ny_mean /= ny[a].size();
      nz_mean /= nz[a].size();
  
      FOR( i, nx[a].size() ) nx_stddev += ((nx[a][i] - nx_mean) * (nx[a][i] - nx_mean));
      FOR( i, ny[a].size() ) ny_stddev += ((ny[a][i] - ny_mean) * (ny[a][i] - ny_mean));
      FOR( i, nz[a].size() ) nz_stddev += ((nz[a][i] - nz_mean) * (nz[a][i] - nz_mean));
  
      nx_stddev = std::sqrt( nx_stddev / (nx[a].size()-1) );
      ny_stddev = std::sqrt( ny_stddev / (ny[a].size()-1) );
      nz_stddev = std::sqrt( nz_stddev / (nz[a].size()-1) );
    }

    mesh.attrF[cnt]     = a;
    mesh.nx_mean[cnt]   = nx_mean;
    mesh.ny_mean[cnt]   = ny_mean;
    mesh.nz_mean[cnt]   = nz_mean;
    mesh.nx_stddev[cnt] = nx_stddev;
    mesh.ny_stddev[cnt] = ny_stddev;
    mesh.nz_stddev[cnt] = nz_stddev;
    cnt++;
  }
  std::cout << "numA_cnt = " << cnt << std::endl;

  msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, centroid, "done" );

  return msu_success;
}

  int 
//---------------------------------------------------------------
  mesh_msu_readBC_triangle
//---------------------------------------------------------------
// read node and element BC information 
(
  Mesh&        mesh, 
  const char*  filebase, 
  int          num_interiorEdgeTags,
  int*         interiorEdgeTags
)
{
  using namespace std;

  assert(filebase);

  bool failed = true;

  std::string bcFilename = std::string(filebase)+".edge";

  std::ifstream bcFile( bcFilename.c_str() );

  if( ! bcFile.good() )
  {
    std::cerr<<"readBC cannot open file " << bcFilename << std::endl;
    failed = true;
    return failed;
  }

  int numF, numA;
  bcFile >> numF >> numA;

  int n, index; std::vector<int> vecN;

  failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN, "init" );

  if( failed ) { std::cerr << "bad init" << std::endl; return failed; }

  // 1. boundary edges

  int indexBase = 1;

  int idx, bc;
  double area, normal[3];

  vecN.resize(2);

  std::set<int> attrF;
  std::map<int,msu_vecD> nx;
  std::map<int,msu_vecD> ny;
  std::map<int,msu_vecD> nz;

  FOR( i, numF )
  {
    // a. set bc 

    bcFile >> idx >> vecN[0] >> vecN[1] >> bc;
    bool itag = 0; FOR(j,num_interiorEdgeTags) if( interiorEdgeTags[j] == bc ) itag = 1;
    if( ! itag )
    {
      for( int j=0; j<3; j++ ) vecN[j] -= indexBase;
      failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN );
      if( failed ) { std::cerr << "did not locate bnd face centroid" << std::endl; return failed; }
      //mesh.elementBoundaryMaterialTypes[index] = flag;
      mesh.elementBoundaryMaterialTypes[index] = bc;
  
      attrF.insert(bc);
  
      // b. compute face normal
  
      msu_mesh_util_area_normal_2d( mesh, vecN, area, normal );
  
      nx[bc].push_back( normal[0] );
      ny[bc].push_back( normal[1] );
      nz[bc].push_back( normal[2] );
    }
  }

  // compute mean and stddev of bc zone face normals

  int nA = attrF.size();
  mesh.num_attrF = nA;
  mesh.attrF     = new int[nA];
  mesh.nx_mean   = new double[nA];
  mesh.ny_mean   = new double[nA];
  mesh.nz_mean   = new double[nA];
  mesh.nx_stddev = new double[nA];
  mesh.ny_stddev = new double[nA];
  mesh.nz_stddev = new double[nA];

  int cnt = 0;
  std::set<int>::iterator I;
  for( I = attrF.begin(); I != attrF.end(); ++I )
  {
    int a = *I;

    double nx_mean=0, nx_stddev=0;
    double ny_mean=0, ny_stddev=0;
    double nz_mean=0, nz_stddev=0;

    FOR( i, nx[a].size() ) nx_mean += nx[a][i];
    FOR( i, ny[a].size() ) ny_mean += ny[a][i];
    FOR( i, nz[a].size() ) nz_mean += nz[a][i];

    nx_mean /= nx[a].size();
    ny_mean /= ny[a].size();
    nz_mean /= nz[a].size();

    FOR( i, nx[a].size() ) nx_stddev += ((nx[a][i] - nx_mean) * (nx[a][i] - nx_mean));
    FOR( i, ny[a].size() ) ny_stddev += ((ny[a][i] - ny_mean) * (ny[a][i] - ny_mean));
    FOR( i, nz[a].size() ) nz_stddev += ((nz[a][i] - nz_mean) * (nz[a][i] - nz_mean));

    nx_stddev = std::sqrt( nx_stddev / (nx[a].size()-1) );
    ny_stddev = std::sqrt( ny_stddev / (ny[a].size()-1) );
    nz_stddev = std::sqrt( nz_stddev / (nz[a].size()-1) );

    mesh.attrF[cnt]     = a;
    mesh.nx_mean[cnt]   = nx_mean;
    mesh.ny_mean[cnt]   = ny_mean;
    mesh.nz_mean[cnt]   = nz_mean;
    mesh.nx_stddev[cnt] = nx_stddev;
    mesh.ny_stddev[cnt] = ny_stddev;
    mesh.nz_stddev[cnt] = nz_stddev;
    cnt++;
  }

  msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN, "done" );

  return msu_success;
}

  int 
//---------------------------------------------------------------
  mesh_msu_readBC_hex
//---------------------------------------------------------------
// read node and element BC information 
(
  Mesh&        mesh, 
  const char*  filebase, 
  int          indexBase
)
{
  using namespace std;

  assert(filebase);

  bool failed = true;

  std::string bcFilename = std::string(filebase)+".bc";

  std::ifstream bcFile( bcFilename.c_str() );

  if( ! bcFile.good() )
  {
    std::cerr<<"readBC cannot open file " << bcFilename << std::endl;
    failed = true;
    return failed;
  }

  std::string firstWord;
  bcFile >> firstWord;
  if( firstWord != "BC" )
  {
    std::cerr<<"firstWord is not BC" << std::endl;
    failed = true;
    return failed;
  }

  int numN, numT, numQ;
  bcFile >> numN >> numT >> numQ;

  int n, index; std::vector<int> vecN;

  failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN, "init" );

  if( failed ) { std::cerr << "bad init" << std::endl; return failed; }

  // 1. boundary quads

  int idx, bc;
  double area, normal[3];

  vecN.resize(4);

  std::set<int> attrF;
  std::map<int,msu_vecD> nx;
  std::map<int,msu_vecD> ny;
  std::map<int,msu_vecD> nz;

  FOR( i, numQ )
  {
    // a. set bc 

    bcFile >> vecN[0] >> vecN[1] >> vecN[2] >> vecN[3] >> bc;
    FOR(j,4)  vecN[j] -= indexBase;
    failed = msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN );
    if( failed ) { std::cerr << "did not locate bnd face centroid" << std::endl; return failed; }
    //mesh.elementBoundaryMaterialTypes[index] = flag;
    mesh.elementBoundaryMaterialTypes[index] = bc;

    attrF.insert(bc);

    // b. compute face normal

    msu_mesh_util_area_normal( mesh, vecN, area, normal );

    nx[bc].push_back( normal[0] );
    ny[bc].push_back( normal[1] );
    nz[bc].push_back( normal[2] );
  }

  // compute mean and stddev of bc zone face normals

  int nA = attrF.size();
  mesh.num_attrF = nA;
  mesh.attrF     = new int[nA];
  mesh.nx_mean   = new double[nA];
  mesh.ny_mean   = new double[nA];
  mesh.nz_mean   = new double[nA];
  mesh.nx_stddev = new double[nA];
  mesh.ny_stddev = new double[nA];
  mesh.nz_stddev = new double[nA];

  int cnt = 0;
  std::set<int>::iterator I;
  for( I = attrF.begin(); I != attrF.end(); ++I )
  {
    int a = *I;

    double nx_mean=0, nx_stddev=0;
    double ny_mean=0, ny_stddev=0;
    double nz_mean=0, nz_stddev=0;

    FOR( i, nx[a].size() ) nx_mean += nx[a][i];
    FOR( i, ny[a].size() ) ny_mean += ny[a][i];
    FOR( i, nz[a].size() ) nz_mean += nz[a][i];

    nx_mean /= nx[a].size();
    ny_mean /= ny[a].size();
    nz_mean /= nz[a].size();

    FOR( i, nx[a].size() ) nx_stddev += ((nx[a][i] - nx_mean) * (nx[a][i] - nx_mean));
    FOR( i, ny[a].size() ) ny_stddev += ((ny[a][i] - ny_mean) * (ny[a][i] - ny_mean));
    FOR( i, nz[a].size() ) nz_stddev += ((nz[a][i] - nz_mean) * (nz[a][i] - nz_mean));

    nx_stddev = std::sqrt( nx_stddev / (nx[a].size()-1) );
    ny_stddev = std::sqrt( ny_stddev / (ny[a].size()-1) );
    nz_stddev = std::sqrt( nz_stddev / (nz[a].size()-1) );

    mesh.attrF[cnt]     = a;
    mesh.nx_mean[cnt]   = nx_mean;
    mesh.ny_mean[cnt]   = ny_mean;
    mesh.nz_mean[cnt]   = nz_mean;
    mesh.nx_stddev[cnt] = nx_stddev;
    mesh.ny_stddev[cnt] = ny_stddev;
    mesh.nz_stddev[cnt] = nz_stddev;
    cnt++;
  }

  msu_mesh_util_find_index_elementBoundaryMaterialTypes( index, mesh, vecN, "done" );

  return msu_success;
}

  int 
//---------------------------------------------------------------
  mesh_msu_read_hex
//---------------------------------------------------------------
// read hex mesh format
( 
  Mesh&       mesh, 
  const char* filebase, 
  int         indexBase
)
{
  using namespace std;
  assert(filebase);
  bool failed=true;
  std::string meshFilename= std::string(filebase)+".mesh";
  std::ifstream meshFile(meshFilename.c_str());

  //std::cout<<"Reading hex mesh: "<<meshFilename<<std::endl;

  if (!meshFile.good())
    {
      std::cerr<<"readHex cannot open file "
         <<meshFilename<<std::endl;
      failed = true;
      return failed;
    }
  //read element type
  std::string fileType;
  meshFile>>fileType;
  if (fileType != "HEX")
    {
      std::cerr<<"readHex does not recognize filetype "
         <<fileType<<std::endl;
      failed = true;
      return failed;
    }
  
  //read mesh size 
  meshFile>>mesh.nNodes_global>>mesh.nElements_global; 

  //read nodes
  mesh.nodeArray         = new double[mesh.nNodes_global*3];
  mesh.nodeMaterialTypes = new int   [mesh.nNodes_global];
  for (int nN=0;nN<mesh.nNodes_global;nN++) 
    {
        meshFile>>mesh.nodeArray[nN*3+0]>>
                  mesh.nodeArray[nN*3+1]>>
                  mesh.nodeArray[nN*3+2];
        
  mesh.nodeMaterialTypes[nN] = 0;    
    }

  //read elements 
  mesh.nNodes_element = 8;
  mesh.elementNodesArray    = new int[mesh.nElements_global*mesh.nNodes_element];
  mesh.elementMaterialTypes = new int[mesh.nElements_global];

  int n0,n1,n2,n3,n4,n5,n6,n7,emt;
  for (int eN=0;eN<mesh.nElements_global;eN++) 
  {
       int eNne = eN*mesh.nNodes_element;
       meshFile>>n0>>n1>>n2>>n3>>n4>>n5>>n6>>n7>>emt;

       mesh.elementNodesArray[eNne+0] = n0-indexBase;
       mesh.elementNodesArray[eNne+1] = n1-indexBase;
       mesh.elementNodesArray[eNne+2] = n2-indexBase;
       mesh.elementNodesArray[eNne+3] = n3-indexBase;
       mesh.elementNodesArray[eNne+4] = n4-indexBase;
       mesh.elementNodesArray[eNne+5] = n5-indexBase;
       mesh.elementNodesArray[eNne+6] = n6-indexBase;
       mesh.elementNodesArray[eNne+7] = n7-indexBase;

       //msu mesh.elementMaterialTypes[eN] = emt-indexBase;
       mesh.elementMaterialTypes[eN] = emt;
  }

  std::string text;
  meshFile >> text >> text >> text;
  meshFile >> mesh.hex_nx >> mesh.hex_ny >> mesh.hex_nz;
 
  return 0;
}
