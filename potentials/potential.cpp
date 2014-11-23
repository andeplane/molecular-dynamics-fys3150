#include <potentials/potential.h>

Potential::Potential() :
    m_potentialEnergy(0)
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

