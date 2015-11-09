#pragma once

class System;
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
    void sample(System *system);
    float sampleKineticEnergy(System *system);
    float samplePotentialEnergy(System *system);
    float sampleTemperature(System *system);
    float sampleDensity(System *system);
    vec3 sampleMomentum(System *system);
};
