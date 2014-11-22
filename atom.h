#pragma once
#include "math/vec3.h"
#include "math/hilbert.h"
#include "math/morton.h"
#include <vector>
using std::vector;

class Atom
{
private:
    static int totalNumberOfAtoms;
    float m_mass;
    unsigned int   m_index;
    unsigned int   m_cellIndex;

public:
    // int   m_numNeighbors;
    // Atom  *m_neighbors[150];
    // vector<Atom*> m_extraNeighbors;

    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom();
    Atom(float mass);
    ~Atom();
    void resetForce();
    void resetVelocityMaxwellian(float temperature);

    float mass() { return m_mass; }
    void setMass(float mass) { m_mass = mass; }
    inline int index() { return m_index; }
    void addNeighbor(Atom *atom);
    void resetNeighbors();
    unsigned int cellIndex() const;
    void setCellIndex(unsigned int cellIndex);
    inline bool operator() (const Atom& atom1, const Atom& atom2)
    {
        // return hilbert_ieee_cmp(3, ) > 0;
        return (atom1.cellIndex() < atom2.cellIndex());
        // return lessThanZOrderDouble((double*)&atom1.position, (double*)&atom2.position);
    }
};
