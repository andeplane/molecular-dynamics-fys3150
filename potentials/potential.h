#pragma once
#include <string>
#include <vector>
#include <system.h>

class Potential
{
protected:
    float m_potentialEnergy;
    float m_pressureVirial;
    unsigned int m_numPairsComputed;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    virtual void calculateForcesAndEnergyAndPressure(System *system) = 0;
    float potentialEnergy();
    float pressureVirial();
    unsigned int numPairsComputed() { return m_numPairsComputed; }
};
