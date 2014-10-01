#pragma once

enum Integrators {VelocityVerlet = 0};

class System;
class Integrator
{
public:
    Integrator();
    virtual ~Integrator() { }
    virtual void integrate(System* system, double timestep) = 0;
};
