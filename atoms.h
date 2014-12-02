#ifndef ATOMS_H
#define ATOMS_H
#include "config.h"

#include "math/vec3.h"


class Atoms
{
public:
    unsigned int numberOfAtoms;
    unsigned int numberOfGhostAtoms;
    MDDataType_t x[MAXNUMATOMS];
    MDDataType_t fx[MAXNUMATOMS];
    MDDataType_t y[MAXNUMATOMS];
    MDDataType_t fy[MAXNUMATOMS];
    MDDataType_t z[MAXNUMATOMS];
    MDDataType_t fz[MAXNUMATOMS];

    int   cellIndex[MAXNUMATOMS];

    MDDataType_t vx[MAXNUMATOMS];
    MDDataType_t vy[MAXNUMATOMS];
    MDDataType_t vz[MAXNUMATOMS];

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
    MDDataType_t x[MAXNUMATOMS];
    MDDataType_t y[MAXNUMATOMS];
    MDDataType_t z[MAXNUMATOMS];
    void update(Atoms &atoms);
    MiniAtoms();
};
#endif // ATOMS_H
