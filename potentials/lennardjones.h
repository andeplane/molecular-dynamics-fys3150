#pragma once
#include <potentials/potential.h>

class LennardJones : public Potential
{
private:
    double m_sigma;
    double m_epsilon;
public:
    LennardJones(double sigma, double epsilon);
    virtual void calculateForces(System *system);
};
