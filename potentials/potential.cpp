#include "potential.h"

Potential::Potential()
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
