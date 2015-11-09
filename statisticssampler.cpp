#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"

StatisticsSampler::StatisticsSampler() :
    m_kineticEnergy(0),
    m_potentialEnergy(0),
    m_temperature(0),
    m_pressure(0)
{

}

void StatisticsSampler::sample(System *system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
}

float StatisticsSampler::sampleKineticEnergy(System *system)
{

}

float StatisticsSampler::samplePotentialEnergy(System *system)
{

}

float StatisticsSampler::sampleTemperature(System *system)
{

}

float StatisticsSampler::sampleDensity(System *system)
{

}

vec3 StatisticsSampler::sampleMomentum(System *system)
{

}
