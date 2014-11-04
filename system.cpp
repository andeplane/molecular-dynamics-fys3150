#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"

System::System() :
    m_potential(0),
    m_integrator(0),
    m_cellList(0),
    m_currentTime(0),
    m_steps(0),
    m_initialized(false)
{
    m_cellList = new CellList();
}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    m_atoms.clear();
}

void System::initialize(float cutoffRadius) {
    m_cellList->setup(this, cutoffRadius);
}

void System::applyPeriodicBoundaryConditions() {
    for(int i=0; i<m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];
        for(int a=0; a<3; a++) {
            if(atom->position[a] < 0) atom->position[a] += m_systemSize[a];
            if(atom->position[a] >= m_systemSize[a]) atom->position[a] -= m_systemSize[a];
        }
    }
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeMomentum() {
    // Initially, when the atoms are given random velocities, there is a non-zero net momentum. We don't want any drift in the system, so we need to remove it.
    StatisticsSampler sampler;
    vec3 momentum = sampler.sampleMomentum(this);
    momentum = momentum / m_atoms.size();
    for(int i=0; i<m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];
        atom->velocity -= momentum;
    }

}

void System::resetForcesOnAllAtoms() {
    for(int i=0; i<m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, float latticeConstant, float temperature) {
    float xCell[4] = {0, 0.5, 0.5, 0};
    float yCell[4] = {0, 0.5, 0, 0.5};
    float zCell[4] = {0, 0, 0.5, 0.5};
    for(int i=0; i< numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j< numberOfUnitCellsEachDimension; j++) {
            for(int k=0; k< numberOfUnitCellsEachDimension; k++) {
                for(int l=0; l<4; l++) {
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    float x = (i+xCell[l])*latticeConstant;
                    float y = (j+yCell[l])*latticeConstant;
                    float z = (k+zCell[l])*latticeConstant;
                    atom->position.set(x,y,z);
                    atom->resetVelocityMaxwellian(temperature);
                    m_atoms.push_back(atom);
                }
            }
        }
    }
    float sideLength = numberOfUnitCellsEachDimension*latticeConstant;
    setSystemSize(vec3(sideLength, sideLength, sideLength));
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(float dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}
