#include <integrators/eulercromer.h>
#include <system.h>

void EulerCromer::integrate(System *system, float timestep)
{
    system->calculateForces();

    system->applyPeriodicBoundaryConditions();
}
