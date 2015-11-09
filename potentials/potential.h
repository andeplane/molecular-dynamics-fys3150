#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <string>
#include <vector>
#include "../system.h"

class Potential
{
protected:
    double m_potentialEnergy = 0;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    double potentialEnergy();
    void setPotentialEnergy(double potentialEnergy);
};
#endif
