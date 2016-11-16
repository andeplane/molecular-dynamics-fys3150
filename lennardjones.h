#ifndef LENNARDJONES_H
#define LENNARDJONES_H

class LennardJones
{
private:
    double m_sigma = 1.0;
    double m_epsilon = 1.0;
    double m_potentialEnergy = 0;

public:
    LennardJones() { }
    void calculateForces(class System &system);
    double potentialEnergy() const;
    double sigma() const;
    void setSigma(double sigma);
    double epsilon() const;
    void setEpsilon(double epsilon);
};
#endif
