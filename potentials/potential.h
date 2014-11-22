#pragma once
#include <string>
#include <vector>
#include <system.h>

class Potential
{
protected:
    float m_potentialEnergy;
    float m_pressureVirial;
    bool m_shouldComputeEnergyAndPressureVirial;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    float potentialEnergy();
    float pressureVirial();

    bool shouldComputeEnergyAndPressureVirial() const;
    void setShouldComputeEnergyAndPressureVirial(bool shouldComputeEnergyAndPressureVirial);
};
