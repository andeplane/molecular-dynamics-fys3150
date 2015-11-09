#pragma once
#include "potential.h"

class LennardJones : public Potential
{
private:
    double m_sigma = 3.405;
    double m_epsilon = 1;
public:
    LennardJones(double sigma, double epsilon);
    ~LennardJones() {}
    virtual void calculateForces(System *system);
};
