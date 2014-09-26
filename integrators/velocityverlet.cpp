#include <integrators/velocityverlet.h>

void VelocityVerlet::firstKick(System *system, double timestep)
{

}

void VelocityVerlet::halfKick(System *system, double timestep)
{

}

void VelocityVerlet::move(System *system, double timestep)
{

}

VelocityVerlet::VelocityVerlet() :
    m_firstStep(false) // This will set the variable m_firstStep to false when the object is created
{

}

void VelocityVerlet::integrate(System *system, double timestep)
{
    if(m_firstStep) {
        firstKick(system, timestep);
        m_firstStep = false;
    }
    halfKick(system, timestep);
    move(system, timestep);
}
