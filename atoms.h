#ifndef ATOMS_H
#define ATOMS_H
#include "config.h"

#include "math/vec3.h"

class MiniAtoms
{
public:
    unsigned int numberOfAtoms;
    unsigned int originalAtom[MAXNUMATOMS];
    float dx[MAXNUMATOMS];
    float dy[MAXNUMATOMS];
    float dz[MAXNUMATOMS];
    float ddx[MAXNUMATOMS];
    float ddy[MAXNUMATOMS];
    float ddz[MAXNUMATOMS];
    float dr2[MAXNUMATOMS];
    float dr2i[MAXNUMATOMS];
    float dr6i[MAXNUMATOMS];

    unsigned long numberOfComputedForces;

    MiniAtoms();
};

class Atoms
{
public:
    unsigned int numberOfAtoms;
    float x[MAXNUMATOMS];
    float fx[MAXNUMATOMS];
    float y[MAXNUMATOMS];
    float fy[MAXNUMATOMS];
    float z[MAXNUMATOMS];
    float fz[MAXNUMATOMS];

    int   cellIndex[MAXNUMATOMS];

    float vx[MAXNUMATOMS];
    float vy[MAXNUMATOMS];
    float vz[MAXNUMATOMS];

    float mass[MAXNUMATOMS];
    int   index[MAXNUMATOMS];

    unsigned long numberOfComputedForces;

    Atoms();
};

#endif // ATOMS_H
