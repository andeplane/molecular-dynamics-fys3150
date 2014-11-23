#pragma once
#include <string>
#include <vector>
#include <system.h>

class Potential
{
protected:
    float m_potentialEnergy;
    float m_pressureVirial;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    virtual void calculateForcesAndEnergyAndPressure(System *system) = 0;
    float potentialEnergy();
    float pressureVirial();
};
