#pragma once
#include "integrators/integrator.h"

class VelocityVerlet : public Integrator
{
public:
    VelocityVerlet() { }
    virtual void integrate(System *system, double dt) override;
};
