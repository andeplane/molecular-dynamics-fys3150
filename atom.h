#pragma once
#include <math/vec3.h>
using CompPhys::vec3;

class Atom
{
private:
    float m_mass;
public:
    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom(float mass);
    ~Atom();
    void resetForce();
    void resetVelocityMaxwellian(float temperature);

    inline float mass() { return m_mass; }
    inline void setMass(float mass) { m_mass = mass; }
};
