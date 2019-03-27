// msu_defines.h

#ifndef msu_defines_h
#define msu_defines_h

#define msu_success 0
#define msu_failure 1

#define FOR(i,n)  for( int i=0; i<(n); ++i )
#define ITER(i,c) for( i=(c).begin(); i!=(c).end(); ++i ) 
#define NINT(a) ((a) >= 0.0 ? (int)((a)+0.5) : (int)((a)-0.5))
#define CONTAINER_HAS_ITEM(C,X) (std::find((C).begin(),(C).end(),X) != (C).end())

#endif
