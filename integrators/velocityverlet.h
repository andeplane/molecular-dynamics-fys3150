#pragma once
#include <integrators/integrator.h>

class VelocityVerlet : public Integrator
{
private:
    void halfKick(System *system, MDDataType_t dt);
    void move(System *system, MDDataType_t dt);
    bool m_firstStep;
public:
    VelocityVerlet();
    virtual void integrate(System *system, MDDataType_t dt);
};
