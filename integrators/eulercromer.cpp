#include <integrators/eulercromer.h>
#include <system.h>

EulerCromer::EulerCromer()
{

}

void EulerCromer::integrate(System *system, double timestep)
{
    system->calculateForces();
    for(Atom *atom : system->atoms()) {
        double timestepDividedByMass = timestep / atom->mass();
        atom->velocity.addAndMultiply(atom->force, timestepDividedByMass); // v += F/m*dt
        atom->position.addAndMultiply(atom->velocity, timestep); // r += v*dt
    }

    system->applyPeriodicBoundaryConditions();
}
