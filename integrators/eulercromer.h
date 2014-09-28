#pragma once
#include <integrators/integrator.h>

class System;
class EulerCromer : public Integrator
{
public:
    EulerCromer();
    virtual void integrate(System* system, double timestep);
};
