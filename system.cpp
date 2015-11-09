#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"

System::System()
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    double xCell[4] = {0, 0.5, 0.5, 0};
    double yCell[4] = {0, 0.5, 0, 0.5};
    double zCell[4] = {0, 0, 0.5, 0.5};
    for(int i=0; i< numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j< numberOfUnitCellsEachDimension; j++) {
            for(int k=0; k< numberOfUnitCellsEachDimension; k++) {
                for(int l=0; l<4; l++) {
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    double x = (i+xCell[l])*latticeConstant;
                    double y = (j+yCell[l])*latticeConstant;
                    double z = (k+zCell[l])*latticeConstant;
                    atom->position.set(x,y,z);
                    atom->resetVelocityMaxwellian(temperature);
                    m_atoms.push_back(atom);
                }
            }
        }
    }
    double sideLength = numberOfUnitCellsEachDimension*latticeConstant;
    setSystemSize(vec3(sideLength, sideLength, sideLength));
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    m_cellList.build(this);
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}
