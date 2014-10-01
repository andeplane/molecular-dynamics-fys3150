#include <integrators/eulercromer.h>
#include <system.h>

void EulerCromer::integrate(System *system, double timestep)
{
    system->calculateForces();
    for(int n=0; n<system->atoms().size(); n++) {
        Atom *atom = system->atoms()[n];
        double timestepDividedByMass = timestep / atom->mass();
        atom->velocity.addAndMultiply(atom->force, timestepDividedByMass); // v += F/m*dt
        atom->position.addAndMultiply(atom->velocity, timestep); // r += v*dt
    }

    system->applyPeriodicBoundaryConditions();
}
