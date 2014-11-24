#include <potentials/potential.h>

Potential::Potential() :
    m_potentialEnergy(0),
    m_pressureVirial(0),
    m_numPairsComputed(0)
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

