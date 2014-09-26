#include <integrators/velocityverlet.h>

VelocityVerlet::VelocityVerlet() :
    m_firstStep(false) // This will set the variable m_firstStep to false when the object is created
{

}

void VelocityVerlet::firstKick(System *system, double dt)
{

}

void VelocityVerlet::halfKick(System *system, double dt)
{

}

void VelocityVerlet::move(System *system, double dt)
{

}

void VelocityVerlet::integrate(System *system, double dt)
{
    if(m_firstStep) {
        firstKick(system, dt);
        m_firstStep = false;
    }
    halfKick(system, dt);
    move(system, dt);

    // Step 1: Fix periodic boundary conditions

    // Step 2: Reset forces on all atoms

    // Step 3: Calculate all forces

    halfKick(system, dt);
}
