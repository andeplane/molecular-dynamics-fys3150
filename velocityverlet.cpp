#include "velocityverlet.h"
#include "system.h"
#include "atom.h"

void VelocityVerlet::integrate(System &system, double dt)
{
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*0.5*dt/atom->mass();
        atom->position += atom->velocity*dt;
    }
    system.applyPeriodicBoundaryConditions();
    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*0.5*dt/atom->mass();
    }
}
