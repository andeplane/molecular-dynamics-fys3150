#pragma once

#ifndef MD_DEBUG
// #define MD_DEBUG
#endif

#ifndef MAXNUMATOMSPERCELL
// Should have length of 4N for simple SIMD
#define MAXNUMATOMSPERCELL 32
#endif

#ifndef CELLSIZE
#define CELLSIZE 2.5
#endif

#ifndef MINIMUMIMAGECONVENTIONTYPE_BRANCH
#define MINIMUMIMAGECONVENTIONTYPE_BRANCH
#endif

#ifndef DISABLEFORCES
// #define DISABLEFORCES
#endif
