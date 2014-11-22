#pragma once
#include <math/vec3.h>
#include <vector>
using std::vector;

class Atom
{
private:
    static int totalNumberOfAtoms;
    float m_mass;
    int   m_index;
public:
    // int   m_numNeighbors;
    // Atom  *m_neighbors[150];
    // vector<Atom*> m_extraNeighbors;

    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom(float mass);
    ~Atom();
    void resetForce();
    void resetVelocityMaxwellian(float temperature);

    float mass() { return m_mass; }
    inline int index() { return m_index; }
    void addNeighbor(Atom *atom);
    void resetNeighbors();
};
