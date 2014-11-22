#include "statisticssampler.h"
#include "potentials/potential.h"
#include "cpelapsedtimer.h"
#include "potentials/lennardjones.h"

StatisticsSampler::StatisticsSampler() :
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_temperature(0),
    m_pressure(0)
{

}

StatisticsSampler::~StatisticsSampler()
{

}

void StatisticsSampler::sample(System *system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    CPElapsedTimer::sampling().start();
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    samplePressure(system);
    CPElapsedTimer::sampling().stop();
}

float StatisticsSampler::sampleKineticEnergy(System *system)
{
    m_kineticEnergy = 0;
    Atoms &atoms = system->atoms();
    for(unsigned int i=0; i<atoms.numberOfAtoms; i++) {
        m_kineticEnergy += 0.5*atoms.mass[i]*(atoms.vx[i]*atoms.vx[i] + atoms.vy[i]*atoms.vy[i] + atoms.vz[i]*atoms.vz[i]);
    }

    return m_kineticEnergy;
}

float StatisticsSampler::samplePotentialEnergy(System *system)
{
    m_potentialEnergy = system->potential()->potentialEnergy();
    return m_potentialEnergy;
}

float StatisticsSampler::sampleTemperature(System *system)
{
    m_temperature = 2.0*m_kineticEnergy/(3*system->atoms().numberOfAtoms);
    return m_temperature;
}

float StatisticsSampler::sampleDensity(System *system)
{
    m_density = system->atoms().numberOfAtoms / system->volume();
    return m_density;
}

float StatisticsSampler::samplePressure(System *system)
{
    float idealGasPressure = m_density*m_temperature;
    float virialPressure = ((LennardJones*)system->potential())->pressureVirial();
//    float virialPressure = 0;
//    for(int i=0; i<system->atoms().size(); i++) {
//        Atom *atom = system->atoms()[i];
//        virialPressure += atom->position.dot(atom->force);
//    }
    virialPressure /= 3*system->volume();
    m_pressure = idealGasPressure + virialPressure;
    return m_pressure;
}

float StatisticsSampler::totalEnergy()
{
    return m_potentialEnergy + m_kineticEnergy;
}

vec3 StatisticsSampler::sampleMomentum(System *system)
{
    m_momentum.setToZero();
    Atoms &atoms = system->atoms();
    for(int i=0; i<system->atoms().numberOfAtoms; i++) {
        m_momentum[0] += atoms.vx[i]*atoms.mass[i];
        m_momentum[1] += atoms.vy[i]*atoms.mass[i];
        m_momentum[2] += atoms.vz[i]*atoms.mass[i];
    }
    return m_momentum;
}
