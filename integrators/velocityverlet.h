#pragma once
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
private:
    void halfKick(System *system, float dt);
    void move(System *system, float dt);
    bool m_firstStep;
public:
    VelocityVerlet();
    virtual void integrate(System *system, float dt);
};
