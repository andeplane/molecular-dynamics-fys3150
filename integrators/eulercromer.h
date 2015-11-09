#ifndef EULERCROMER_H
#define EULERCROMER_H
#include "integrators/integrator.h"

class System;
class EulerCromer : public Integrator
{
public:
    EulerCromer() {}
    ~EulerCromer() {}
    virtual void integrate(System* system, double timestep);
};

#endif
