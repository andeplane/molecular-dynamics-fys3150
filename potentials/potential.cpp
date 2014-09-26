#include <potentials/potential.h>

Potential::Potential() :
    m_potentialEnergy(0)
{

}

double Potential::potentialEnergy()
{
    return m_potentialEnergy;
}

void Potential::setPotentialEnergy(double potentialEnergy)
{
    m_potentialEnergy = potentialEnergy;
}
