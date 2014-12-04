#pragma once
#include "potentials/potential.h"

class LennardJones : public Potential
{
private:
    MDDataType_t m_sigma;
    MDDataType_t m_sigma6;
    MDDataType_t m_epsilon;
    MDDataType_t m_24epsilon;
    MDDataType_t m_rCutSquared;
    MDDataType_t m_potentialEnergyAtRcut;
    void calculateForcesAllPairs(System *system);
public:
    LennardJones(MDDataType_t sigma, MDDataType_t epsilon, MDDataType_t cutoffRadius);
    ~LennardJones() {}
    virtual void calculateForces(System *system);
    virtual void calculateForcesAndEnergyAndPressure(System *system);
    void calculateForcesBetweenAtoms(Atom *atom1, Atom *atom2, vec3 &deltaRVector, const vec3 &systemSize);
};
