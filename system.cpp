#include <system.h>
#include <integrators/integrator.h>
#include <potentials/potential.h>

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0)
{

}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::resetForcesOnAllAtoms() {

}

void System::createFCCLattice(int numberOfUnitCells, double latticeConstant)
{

}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}
