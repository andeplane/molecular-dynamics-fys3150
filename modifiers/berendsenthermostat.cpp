#include "berendsenthermostat.h"
#include "statisticssampler.h"
#include "system.h"
#include "atom.h"
#include <cmath>

BerendsenThermostat::BerendsenThermostat(float temperature, float relaxationFactor) :
    m_temperature(temperature),
    m_relaxationFactor(relaxationFactor)
{
}

void BerendsenThermostat::apply(System *system, StatisticsSampler *sampler)
{
    float scaleFactor = sqrt(1 + m_relaxationFactor*(sampler->temperature()/m_temperature - 1));
    for(Atom &atom : system->atoms()) {
        atom.velocity *= scaleFactor;
    }
}
