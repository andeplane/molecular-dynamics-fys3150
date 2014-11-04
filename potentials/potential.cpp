#include <potentials/potential.h>

Potential::Potential() :
    m_potentialEnergy(0)
{

}

float Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

void Potential::setPotentialEnergy(float potentialEnergy)
{
    m_potentialEnergy = potentialEnergy;
}

void Potential::addPotentialEnergy(float potentialEnergy)
{
    m_potentialEnergy = +potentialEnergy;
}
