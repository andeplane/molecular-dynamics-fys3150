#ifndef ATOMS_H
#define ATOMS_H
#define MAXNUMATOMS 100000
class Atoms
{
public:
    Atoms();
    unsigned int numberOfAtoms;
    float x[MAXNUMATOMS];
    float y[MAXNUMATOMS];
    float z[MAXNUMATOMS];

    float fx[MAXNUMATOMS];
    float fy[MAXNUMATOMS];
    float fz[MAXNUMATOMS];

    float vx[MAXNUMATOMS];
    float vy[MAXNUMATOMS];
    float vz[MAXNUMATOMS];

    float mass[MAXNUMATOMS];
    int   index[MAXNUMATOMS];
};

#endif // ATOMS_H
