#pragma once
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
private:
    void firstKick(System *system, double timestep);
    void halfKick(System *system, double timestep);
    void move(System *system, double timestep);
    bool m_firstStep;
public:
    VelocityVerlet();

    virtual void integrate(System *system, double timestep);
};
