#include <cmath>
#include <atom.h>
#include <math/random.h>

Atom::Atom(double mass) :
    m_mass(mass)
{
    
}

Atom::~Atom()
{

}

void Atom::resetForce()
{
    force.setToZero();
}

void Atom::resetVelocityMaxwellian(double temperature)
{
    // Resetting the velocity according to a Maxwell-Boltzmann distribution (see http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution )
    double boltzmannConstant = 1.0; // In atomic units, the boltzmann constant equals 1
    double standardDeviation = sqrt(boltzmannConstant*temperature/m_mass);
    velocity.randomGaussian(0, standardDeviation);
}
