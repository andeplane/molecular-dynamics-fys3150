#include "atoms.h"
#include "math/hilbert.h"
#include <algorithm>
#include <cstring>

#include <iostream>

Atoms::Atoms() :
    numberOfAtoms(0),
    numberOfComputedForces(0)
{
    memset(cellIndex,0,MAXNUMATOMS*sizeof(int));
    memset(index,0,MAXNUMATOMS*sizeof(int));
}

MiniAtoms::MiniAtoms() :
    numberOfAtoms(0),
    numberOfComputedForces(0)
{

}
