#include "atoms.h"
#include "math/hilbert.h"
#include <algorithm>
#include <cstring>

#include <iostream>

Atoms::Atoms() :
    numberOfAtoms(0),
    numberOfGhostAtoms(0),
    numberOfComputedForces(0)
{
    memset(cellIndex,0,MAXNUMATOMS*sizeof(int));
    memset(index,0,MAXNUMATOMS*sizeof(int));
}

void Atoms::sort()
{
    for (int i = 1; i < numberOfAtoms; i++) {
        int j = i;

        while (j > 0 && cellIndex[j-1] > cellIndex[j]) {
            std::swap(x[j], x[j-1]);
            std::swap(y[j], y[j-1]);
            std::swap(z[j], z[j-1]);

            std::swap(vx[j], vx[j-1]);
            std::swap(vy[j], vy[j-1]);
            std::swap(vz[j], vz[j-1]);

            std::swap(fx[j], fx[j-1]);
            std::swap(fy[j], fy[j-1]);
            std::swap(fz[j], fz[j-1]);

            std::swap(mass[j], mass[j-1]);
            std::swap(index[j], index[j-1]);
            std::swap(cellIndex[j], cellIndex[j-1]);

            j--;
        }//end of while loop
    }//end of for loop
}
