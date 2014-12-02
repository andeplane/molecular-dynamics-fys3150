#pragma once
#include "config.h"

class System;
class Integrator
{
public:
    Integrator();
    virtual ~Integrator() { }
    virtual void integrate(System* system, MDDataType_t timestep) = 0;
};
