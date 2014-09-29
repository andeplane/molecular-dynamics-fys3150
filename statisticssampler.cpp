#include <statisticssampler.h>

StatisticsSampler::StatisticsSampler()
{

}

StatisticsSampler::~StatisticsSampler()
{

}

void StatisticsSampler::sample(System *system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    // ...
}

void StatisticsSampler::sampleKineticEnergy(System *system)
{

}
