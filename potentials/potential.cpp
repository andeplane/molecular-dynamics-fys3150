#include <potentials/potential.h>

Potential::Potential() :
    m_potentialEnergy(0)
{
    
}

MDDataType_t Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

MDDataType_t Potential::pressureVirial()
{
    return m_pressureVirial;
}

