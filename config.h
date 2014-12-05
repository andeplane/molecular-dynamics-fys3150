#pragma once

#ifndef SINGLEPRECISION
// #define SINGLEPRECISION
#endif

#ifdef SINGLEPRECISION
typedef float MDDataType_t;
#else
typedef double MDDataType_t;
#endif

#ifndef CORRECTNEIGHBORLISTUPDATE
#define CORRECTNEIGHBORLISTUPDATE
#endif

#ifndef MAXNUMATOMS
#define MAXNUMATOMS 200000
#endif

#ifndef MAXNUMNEIGHBORS
#define MAXNUMNEIGHBORS 300
#endif

#ifndef BENCHMARK
// #define BENCHMARK
#endif

#ifndef BUILDNEIGHBORLIST
#define BUILDNEIGHBORLIST 20
#endif
