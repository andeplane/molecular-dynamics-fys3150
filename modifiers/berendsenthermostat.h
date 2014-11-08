#ifndef BERENDSENTHERMOSTAT_H
#define BERENDSENTHERMOSTAT_H
class System;
class StatisticsSampler;

class BerendsenThermostat
{
private:
    float m_temperature;
    float m_relaxationFactor;
public:
    BerendsenThermostat(float temperature, float relaxationFactor);
    void apply(System *system, StatisticsSampler *sampler);
};

#endif // BERENDSENTHERMOSTAT_H
