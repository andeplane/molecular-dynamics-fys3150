#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {

    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{

}

void StatisticsSampler::sampleTemperature(System &system)
{

}

void StatisticsSampler::sampleDensity(System &system)
{

}
