#ifndef ATOMS_H
#define ATOMS_H
#include "config.h"

#include "math/vec3.h"


class Atoms
{
public:
    unsigned int numberOfAtoms;
    unsigned int numberOfGhostAtoms;
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

    float inverseMass[MAXNUMATOMS];
    int   index[MAXNUMATOMS];

    unsigned long numberOfComputedForces;

    Atoms();
    void sort();
    unsigned int numberOfAtomsIncludingGhosts() { return numberOfAtoms + numberOfGhostAtoms; }
};

class MiniAtoms
{
public:
    unsigned int numberOfAtoms;
    float x[MAXNUMATOMS];
    float y[MAXNUMATOMS];
    float z[MAXNUMATOMS];
    void update(Atoms &atoms);
    MiniAtoms();
};
#endif // ATOMS_H
