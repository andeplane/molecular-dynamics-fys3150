#pragma once
#include <math/vec3.h>
using CompPhys::vec3;

class Atom
{
private:
    double m_mass;
public:
    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom(double mass);
    ~Atom();
    void resetForce();
    void resetVelocityMaxwellian(double temperature);

    inline double mass() { return m_mass; }
    inline void setMass(double mass) { m_mass = mass; }
};
