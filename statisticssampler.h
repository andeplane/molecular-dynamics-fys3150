#pragma once
#include <system.h>

class StatisticsSampler
{
private:
    vec3  m_momentum;
    float m_kineticEnergy;
    float m_potentialEnergy;
    float m_temperature;
    float m_pressure;
    float m_density;
public:
    StatisticsSampler();
    ~StatisticsSampler();
    void  sample(System *system);
    float sampleKineticEnergy(System *system);
    float samplePotentialEnergy(System *system);
    float sampleTemperature(System *system);
    float sampleDensity(System *system);
    float samplePressure(System *system);
    vec3 sampleMomentum(System *system);

    vec3 momentum() { return m_momentum; }
    float totalEnergy();
    float kineticEnergy() { return m_kineticEnergy; }
    float potentialEnergy() { return m_potentialEnergy; }
    float temperature() { return m_temperature; }
    float pressure() { return m_pressure; }
    float density() { return m_density; }
};
