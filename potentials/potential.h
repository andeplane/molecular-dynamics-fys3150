#pragma once
#include <string>
#include <vector>
#include <system.h>

class Potential
{
protected:
    MDDataType_t m_potentialEnergy;
    MDDataType_t m_pressureVirial;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    virtual void calculateForcesAndEnergyAndPressure(System *system) = 0;
    MDDataType_t potentialEnergy();
    MDDataType_t pressureVirial();
};
