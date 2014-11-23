#ifndef ATOMS_H
#define ATOMS_H
#define MAXNUMATOMS 50000
#include "math/vec3.h"

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
    void sort();
};

#endif // ATOMS_H
