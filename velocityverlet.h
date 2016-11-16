#ifndef VELOCITYVERLET_H
#define VELOCITYVERLET_H

class VelocityVerlet
{
private:
    bool m_firstStep = true;
public:
    VelocityVerlet() {}
    void integrate(class System &system, double dt);
};
#endif
