#include "atoms.h"
#include "math/morton.h"
#include "math/hilbert.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <iostream>
#include <numeric> // iota
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
    return;

    std::vector<unsigned int> idx(numberOfAtoms);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int i, int j) {
        double a[3] = { x[i], y[i], z[i] };
        double b[3] = { x[j], y[j], z[j] };
#ifdef MD_HILBERT
        return hilbert_ieee_cmp(3, a, b) < 0;
#else
        return lessThanZOrderDouble(a, b);
#endif
    });

#define COPY(var) float old_##var[MAXNUMATOMS]; memcpy(old_##var, var, numberOfAtoms*sizeof(float));
    COPY(x); COPY(y); COPY(z);
    COPY(vx); COPY(vy); COPY(vz);
#undef COPY

    std::transform(idx.begin(), idx.end(), x, [&](int i) { return old_x[i]; });
    std::transform(idx.begin(), idx.end(), y, [&](int i) { return old_y[i]; });
    std::transform(idx.begin(), idx.end(), z, [&](int i) { return old_z[i]; });
    std::transform(idx.begin(), idx.end(), vz, [&](int i) { return old_vx[i]; });
    std::transform(idx.begin(), idx.end(), vy, [&](int i) { return old_vy[i]; });
    std::transform(idx.begin(), idx.end(), vz, [&](int i) { return old_vz[i]; });
}


void MiniAtoms::update(Atoms &atoms)
{
    memcpy(x,atoms.x,sizeof(x));
    memcpy(y,atoms.y,sizeof(y));
    memcpy(z,atoms.z,sizeof(z));
}

MiniAtoms::MiniAtoms() :
    numberOfAtoms(0)
{

}
