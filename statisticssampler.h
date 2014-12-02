#pragma once
#include "system.h"
#include "io.h"
class StatisticsSampler
{
private:
    IO   *m_fileHandler;
    vec3  m_momentum;
    double m_kineticEnergy;
    double m_potentialEnergy;
    double m_temperature;
    double m_pressure;
    double m_density;
public:
    StatisticsSampler(IO *fileHandler);
    ~StatisticsSampler();
    void  sample(System *system);
    double sampleKineticEnergy(System *system);
    double samplePotentialEnergy(System *system);
    double sampleTemperature(System *system);
    double sampleDensity(System *system);
    double samplePressure(System *system);
    vec3 sampleMomentum(System *system);

    vec3 momentum() { return m_momentum; }
    double totalEnergy();
    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double temperature() { return m_temperature; }
    double pressure() { return m_pressure; }
    double density() { return m_density; }
};
