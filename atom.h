#pragma once
#include <math/vec3.h>

class Atom
{
private:
    static int totalNumberOfAtoms;
    float m_mass;
    int   m_index;
public:
    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom(float mass);
    ~Atom();
    void resetForce();
    void resetVelocityMaxwellian(float temperature);

    float mass() { return m_mass; }
    inline int index() { return m_index; }
};
