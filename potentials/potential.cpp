#include <potentials/potential.h>


bool Potential::shouldComputeEnergyAndPressureVirial() const
{
    return m_shouldComputeEnergyAndPressureVirial;
}

void Potential::setShouldComputeEnergyAndPressureVirial(bool shouldComputeEnergyAndPressureVirial)
{
    m_shouldComputeEnergyAndPressureVirial = shouldComputeEnergyAndPressureVirial;
}
Potential::Potential() :
    m_potentialEnergy(0),
    m_shouldComputeEnergyAndPressureVirial(false)
{
    
}

float Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

float Potential::pressureVirial()
{
    return m_pressureVirial;
}

