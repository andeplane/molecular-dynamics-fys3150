#include "statisticssampler.h"
#include "potentials/potential.h"
#include "cpelapsedtimer.h"
#include "potentials/lennardjones.h"

StatisticsSampler::StatisticsSampler(IO *fileHandler) :
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_temperature(0),
    m_pressure(0),
    m_fileHandler(fileHandler)
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

    if(m_fileHandler) {
        m_fileHandler->writeStatistics(system->currentTime(), m_kineticEnergy, m_potentialEnergy, m_pressure, m_temperature);
    }
    CPElapsedTimer::sampling().stop();
}

double StatisticsSampler::sampleKineticEnergy(System *system)
{
    m_kineticEnergy = 0;
    Atoms &atoms = system->atoms();
    // double kineticEnergy = 0;
//#ifdef MD_SIMD
//#pragma simd reduction(+: kineticEnergy)
//#endif
    for(unsigned int i=0; i<atoms.numberOfAtoms; i++) {
        m_kineticEnergy += 0.5/atoms.inverseMass[i]*(atoms.vx[i]*atoms.vx[i] + atoms.vy[i]*atoms.vy[i] + atoms.vz[i]*atoms.vz[i]);
    }

    // m_kineticEnergy = kineticEnergy;
    return m_kineticEnergy;
}

double StatisticsSampler::samplePotentialEnergy(System *system)
{
    m_potentialEnergy = system->potential()->potentialEnergy();
    return m_potentialEnergy;
}

double StatisticsSampler::sampleTemperature(System *system)
{
    m_temperature = 2.0*m_kineticEnergy/(3*system->atoms().numberOfAtoms);
    return m_temperature;
}

double StatisticsSampler::sampleDensity(System *system)
{
    m_density = system->atoms().numberOfAtoms / system->volume();
    return m_density;
}

double StatisticsSampler::samplePressure(System *system)
{
    double idealGasPressure = m_density*m_temperature;
    double virialPressure = ((LennardJones*)system->potential())->pressureVirial();
//    double virialPressure = 0;
//    for(int i=0; i<system->atoms().size(); i++) {
//        Atom *atom = system->atoms()[i];
//        virialPressure += atom->position.dot(atom->force);
//    }
    virialPressure /= 3*system->volume();
    m_pressure = idealGasPressure + virialPressure;
    return m_pressure;
}

double StatisticsSampler::totalEnergy()
{
    return m_potentialEnergy + m_kineticEnergy;
}

vec3 StatisticsSampler::sampleMomentum(System *system)
{
    m_momentum.setToZero();
    Atoms &atoms = system->atoms();
    for(int i=0; i<system->atoms().numberOfAtoms; i++) {
        m_momentum[0] += atoms.vx[i]/atoms.inverseMass[i];
        m_momentum[1] += atoms.vy[i]/atoms.inverseMass[i];
        m_momentum[2] += atoms.vz[i]/atoms.inverseMass[i];
    }
    return m_momentum;
}
