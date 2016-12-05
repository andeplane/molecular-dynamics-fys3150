#include "statisticssampler.h"
#include "potentials/potential.h"
#include "cpelapsedtimer.h"
#include "potentials/lennardjones.h"
#include <iostream>
#include <iomanip>
using namespace std;

string StatisticsSampler::filename() const
{
    return m_filename;
}

void StatisticsSampler::setFilename(const string &filename)
{
    m_filename = filename;
}

StatisticsSampler::StatisticsSampler() :
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_temperature(0),
    m_pressure(0),
    m_density(0),
    m_msd(0)
{
    setFilename("statistics.txt");
}

StatisticsSampler::~StatisticsSampler()
{
    m_file.close();
}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    if(!m_file.good()) {
        m_file.open(m_filename, std::ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            std::cout << "Error, could not open " <<  m_filename << std::endl;
            exit(1);
        }
        cout << "Opened file " << m_filename << endl;
        m_file << setw(20) << "Timestep" <<
                    setw(20) << "Time" <<
                    setw(20) << "Temperature" <<
                    setw(20) << "KineticEnergy" <<
                    setw(20) << "PotentialEnergy" <<
                    setw(20) << "TotalEnergy" <<
                    setw(20) << "Pressure" <<
                    setw(20) << "MSD" << endl;
    }

    m_file << system.steps() <<
        setw(20) << system.currentTime() <<
        setw(20) << temperature() <<
        setw(20) << kineticEnergy() <<
        setw(20) << potentialEnergy() <<
        setw(20) << totalEnergy() <<
        setw(20) << pressure() <<
        setw(20) << msd() << endl;
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
    sampleMSD(system);
    saveToFile(*system);
    CPElapsedTimer::sampling().stop();
}

double StatisticsSampler::sampleKineticEnergy(System *system)
{
    m_kineticEnergy = 0;
    for(Atom &atom : system->atoms()) {
        m_kineticEnergy += 0.5*atom.mass()*atom.velocity.lengthSquared();
    }

    return m_kineticEnergy;
}

double StatisticsSampler::samplePotentialEnergy(System *system)
{
    m_potentialEnergy = system->potential()->potentialEnergy();
    return m_potentialEnergy;
}

double StatisticsSampler::sampleTemperature(System *system)
{
    m_temperature = 2.0*m_kineticEnergy/(3*system->atoms().size());
    return m_temperature;
}

double StatisticsSampler::sampleDensity(System *system)
{
    m_density = system->atoms().size() / system->volume();
    return m_density;
}

double StatisticsSampler::samplePressure(System *system)
{
    float idealGasPressure = m_density*m_temperature;
    float virialPressure = ((LennardJones*)system->potential())->pressureVirial();
    virialPressure /= 3*system->volume();
    m_pressure = idealGasPressure + virialPressure;
    return m_pressure;
}

double StatisticsSampler::totalEnergy()
{
    return m_potentialEnergy + m_kineticEnergy;
}

void StatisticsSampler::sampleMSD(System *system)
{
    m_msd = 0;
    for(Atom &atom : system->atoms()) {
        m_msd += (atom.position - atom.initialPosition).lengthSquared();
    }
    m_msd /= system->atoms().size();
}

vec3 StatisticsSampler::sampleMomentum(System *system)
{
    m_momentum.setToZero();

    for(Atom &atom : system->atoms()) {
        m_momentum.addAndMultiply(atom.velocity, atom.mass());
    }
    return m_momentum;
}
