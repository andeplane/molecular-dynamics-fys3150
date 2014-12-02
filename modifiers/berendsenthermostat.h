#ifndef BERENDSENTHERMOSTAT_H
#define BERENDSENTHERMOSTAT_H
class System;
class StatisticsSampler;

class BerendsenThermostat
{
private:
    double m_temperature;
    double m_relaxationFactor;
public:
    BerendsenThermostat(double temperature, double relaxationFactor);
    void apply(System *system, StatisticsSampler *sampler);
};

#endif // BERENDSENTHERMOSTAT_H
