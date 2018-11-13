// msu_mesh_util.ipp

using namespace std;

  static void 
//-----------------------------------------------------------
  s_cross_product
//-----------------------------------------------------------
(
  double* n1,
  double* n2,
  double* n3,
  double* prod
)
{
  double x1,x2,x3, y1,y2,y3, z1,z2,z3;
  double x21, y21, z21;
  double x31, y31, z31;

  x1 = n1[0]; y1 = n1[1]; z1 = n1[2];
  x2 = n2[0]; y2 = n2[1]; z2 = n2[2];
  x3 = n3[0]; y3 = n3[1]; z3 = n3[2];
 
  x21 = x2-x1; x31 = x3-x1;
  y21 = y2-y1; y31 = y3-y1;
  z21 = z2-z1; z31 = z3-z1;

  prod[0] =  y21*z31 - z21*y31;
  prod[1] = -x21*z31 + z21*x31;
  prod[2] =  x21*y31 - y21*x31;
}

  static void 
//-----------------------------------------------------------
  s_face_area_normal
//-----------------------------------------------------------
// area always positive
(
  const msu_vecI&  vecN,
  double*          nodeArray,
  double&          area,
  double*          normal 
)
{
  int numN = vecN.size();

  double zero = 1.e-30;

  area = 0.;
  FOR(i,3) normal[i] = 0.;

  if( numN == 3 )
  {
    int i0 = vecN[0]*3;
    int i1 = vecN[1]*3;
    int i2 = vecN[2]*3;
    s_cross_product( &nodeArray[i0], &nodeArray[i1], &nodeArray[i2], normal );
    double mag = std::sqrt( normal[0]*normal[0] +
                            normal[1]*normal[1] +
                            normal[2]*normal[2] );
    FOR(i,3) normal[i] /= (mag + zero);
    area = 0.5 * mag; // always (+)
  }
  else
  {
    // no guarentee that the face normal is correct.

    double fn[3];
    double fc[3]; 
    FOR(i,3) fc[i] = 0;
    FOR( n, numN )
    {
      FOR(i,3) fc[i] += nodeArray[vecN[n]*3+i]; 
    } 
    FOR(i,3) fc[i] /= numN;
  
    FOR( i, numN )
    {
      int j = i + 1;
      if( j == numN ) j = 0;
      double fa;
      int i0 = vecN[i]*3;
      int i1 = vecN[j]*3;
      s_cross_product( &nodeArray[i0], &nodeArray[i1], fc, fn );
      double mag = std::sqrt( fn[0]*fn[0] +
                              fn[1]*fn[1] +
                              fn[2]*fn[2] );
      FOR(j,3) fn[j] /= (mag + zero);
      fa = 0.5 * mag; // always (+)
      area += fa;
      FOR(j,3) normal[j] += fn[j]*fa;
    }
    FOR(i,3) normal[i] /= area; // area-weighted normal, always (+)
  }
}

  static void 
//-----------------------------------------------------------
  s_face_area_normal
//-----------------------------------------------------------
// area always positive
(
  const msu_vecI&  vecN,
  msu_vecD&        vecX,
  msu_vecD&        vecY,
  msu_vecD&        vecZ,
  double&          area,
  double*          normal 
)
{
  int numN = vecN.size();

  double zero = 1.e-30;

  area = 0.;
  FOR(i,3) normal[i] = 0.;

  vector<msu_vecD> pos(numN);

  FOR( i, numN ) 
  {
    pos[i].resize(3);
    pos[i][0] = vecX[vecN[i]];
    pos[i][1] = vecY[vecN[i]];
    pos[i][2] = vecZ[vecN[i]];
  }

  if( numN == 3 )
  {
    s_cross_product( &pos[0][0], &pos[1][0], &pos[2][0], normal );
    double mag = std::sqrt( normal[0]*normal[0] +
                            normal[1]*normal[1] +
                            normal[2]*normal[2] );
    FOR(i,3) normal[i] /= (mag + zero);
    area = 0.5 * mag; // always (+)
  }
  else
  {
    // no guarentee that the face normal is correct.

    double fn[3];
    double fc[3]; 
    FOR(i,3) fc[i] = 0;
    FOR( n, numN )
    {
      FOR(i,3) fc[i] += pos[vecN[n]][i]; 
    } 
    FOR(i,3) fc[i] /= numN;
  
    FOR( i, numN )
    {
      int j = i + 1;
      if( j == numN ) j = 0;
      double fa;
      s_cross_product( &pos[0][0], &pos[1][0], fc, fn );
      double mag = std::sqrt( fn[0]*fn[0] +
                              fn[1]*fn[1] +
                              fn[2]*fn[2] );
      FOR(j,3) fn[j] /= (mag + zero);
      fa = 0.5 * mag; // always (+)
      area += fa;
      FOR(j,3) normal[j] += fn[j]*fa;
    }
    FOR(i,3) normal[i] /= area; // area-weighted normal, always (+)
  }
}