#pragma once
#include "potentials/potential.h"

class LennardJones : public Potential
{
private:
    float m_sigma;
    float m_sigma6;
    float m_epsilon;
    float m_rCutSquared;
    float m_potentialEnergyAtRcut;
    int   m_timeSinceLastNeighborListUpdate;
public:
    LennardJones(float sigma, float epsilon, float cutoffRadius);
    ~LennardJones() {}
    virtual void calculateForces(System *system);
    void calculateForcesBetweenAtoms(Atom *atom1, Atom *atom2, vec3 &deltaRVector, const vec3 &systemSize);
};
