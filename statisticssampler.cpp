#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>

using std::ofstream; using std::cout; using std::endl;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system)
{
    // Save the statistical properties for each timestep for plotting etc.
    // First, open the file if it's not open already
    if(!m_file.good()) {
        m_file.open("statistics.txt", ofstream::out);
        // If it's still not open, something bad happened...
        if(!m_file.good()) {
            cout << "Error, could not open statistics.txt" << endl;
            exit(1);
        }
    }

    // Print out values here
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
    m_potentialEnergy = system.potential().potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    // Hint: reuse the kinetic energy that we already calculated
}

void StatisticsSampler::sampleDensity(System &system)
{

}
