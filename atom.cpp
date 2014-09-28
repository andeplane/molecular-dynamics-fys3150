#include <cmath>
#include <atom.h>
#include <math/random.h>

Atom::Atom(double mass) :
    m_mass(mass)
{
    
}

void Atom::resetForce()
{
    force.setToZero();
}

void Atom::resetVelocityMaxwellian(double temperature)
{
    // Resetting the velocity according to a Maxwell-Boltzmann distribution (see http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution )
    double standardDeviation = sqrt(temperature/m_mass);
    velocity.x = Random::nextGaussian(0, standardDeviation);
    velocity.y = Random::nextGaussian(0, standardDeviation);
    velocity.z = Random::nextGaussian(0, standardDeviation);
}
