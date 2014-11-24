#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "cpelapsedtimer.h"
#include <stdlib.h>
#include <cmath>
#include <cassert>

using namespace std;

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0),
    m_initialized(false),
    m_shouldSample(false),
    numberOfAtoms(0)
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
}

void System::validate() {
    if(m_cellList.numberOfCellsX() < 3 || m_cellList.numberOfCellsY() < 3 || m_cellList.numberOfCellsZ() < 3) {
        cout << "Error, system size too small to have at least 3 cells in each dimension, aborting." << endl;
        cout << "Minimum system size: " << UnitConverter::lengthToAngstroms(3*m_rShell) << " Å" << endl;
        exit(1);
    }
}

void System::printStatus() {
    cout << "System size: " << m_systemSize << endl;
    cout << "Number of cells: " << m_cellList.numberOfCellsX() << ", " << m_cellList.numberOfCellsY() << ", " << m_cellList.numberOfCellsZ() << endl;
    cout << "Cell size: " << UnitConverter::lengthToAngstroms(m_cellList.lengthX()) << " x " << UnitConverter::lengthToAngstroms(m_cellList.lengthY()) << " x " << UnitConverter::lengthToAngstroms(m_cellList.lengthZ()) << " Å^3" << endl;
}

void System::applyPeriodicBoundaryConditions() {
    CPElapsedTimer::periodicBoundaryConditions().start();

    m_cellList.forEachAtom([&](Cell &cell, unsigned int i) {
        // Fixes the bug where the position is slightly negative, but due to float precision, both tests (<0) and >=size) pass
        if(cell.x[i] < 0 && cell.x[i] > -1e-6) cell.x[i] = 0;
        if(cell.y[i] < 0 && cell.y[i] > -1e-6) cell.y[i] = 0;
        if(cell.z[i] < 0 && cell.z[i] > -1e-6) cell.z[i] = 0;

        if(cell.x[i] < 0) cell.x[i] += m_systemSize[0];
        if(cell.x[i] >= m_systemSize[0]) cell.x[i] -= m_systemSize[0];
        if(cell.y[i] < 0) cell.y[i] += m_systemSize[1];
        if(cell.y[i] >= m_systemSize[1]) cell.y[i] -= m_systemSize[1];
        if(cell.z[i] < 0) cell.z[i] += m_systemSize[2];
        if(cell.z[i] >= m_systemSize[2]) cell.z[i] -= m_systemSize[2];
    });

    CPElapsedTimer::periodicBoundaryConditions().stop();
    m_cellList.update();
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
}

void System::removeMomentum() {
    // Initially, when the atoms are given random velocities, there is a non-zero net momentum. We don't want any drift in the system, so we need to remove it.
    StatisticsSampler sampler;
    vec3 momentum = sampler.sampleMomentum(this);
    
    momentum /= numberOfAtoms;
    m_cellList.forEachAtom([&](Cell &cell, unsigned int i) {
        cell.vx[i] -= momentum[0]/cell.mass[i];
        cell.vy[i] -= momentum[1]/cell.mass[i];
        cell.vz[i] -= momentum[2]/cell.mass[i];
    });
}

void System::resetForcesOnAllAtoms() {
    m_cellList.forEachAtom([&](Cell &cell, unsigned int i) {
        cell.fx[i] = 0;
        cell.fy[i] = 0;
        cell.fz[i] = 0;
    });
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, float latticeConstant, float temperature) {
    float sideLength = numberOfUnitCellsEachDimension*latticeConstant;
    setSystemSize(vec3(sideLength, sideLength, sideLength));

    m_cellList.setup(this);
    printStatus();
    float xCell[4] = {0, 0.5, 0.5, 0};
    float yCell[4] = {0, 0.5, 0, 0.5};
    float zCell[4] = {0, 0, 0.5, 0.5};

    for(int i=0; i< numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j< numberOfUnitCellsEachDimension; j++) {
            for(int k=0; k< numberOfUnitCellsEachDimension; k++) {
                for(int l=0; l<4; l++) {
                    float mass = UnitConverter::massFromSI(6.63352088e-26);
                    float x = (i+xCell[l])*latticeConstant;
                    float y = (j+yCell[l])*latticeConstant;
                    float z = (k+zCell[l])*latticeConstant;

                    vec3 velocity;
                    float boltzmannConstant = 1.0; // In atomic units, the boltzmann constant equals 1
                    float standardDeviation = sqrt(boltzmannConstant*temperature/mass);
                    velocity.randomGaussian(0, standardDeviation);
                    m_cellList.addAtom(vec3(x,y,z), velocity, mass, numberOfAtoms);
                    numberOfAtoms++;
                }
            }
        }
    }

    cout << "System initialized with " << numberOfAtoms << " atoms." << endl;
    m_initialized = true;

    validate();
}

void System::setShouldSample(bool shouldSample)
{
    m_shouldSample = shouldSample;
}

void System::createGhostAtoms()
{
    //    m_atoms->numberOfGhostAtoms = 0; // Reset all ghosts
    
    //    for(unsigned short dimension=0; dimension<3; dimension++) {
    //        unsigned int newGhosts = m_atoms->numberOfGhostAtoms;
    
    //        for(unsigned int i=0; i<m_atoms->numberOfAtoms+newGhosts; i++) {
    //            bool shouldCopy = false;
    //            unsigned int atomIndex = m_atoms->numberOfAtomsIncludingGhosts();
    
    //            if(dimension==0) {
    //                if(m_atoms->x[i] < m_rShell) {
    //                    m_atoms->x[atomIndex] = m_atoms->x[i] + m_systemSize[0]; shouldCopy = true;
    //                } else if (m_atoms->x[i] > m_systemSize[0] - m_rShell) {
    //                    m_atoms->x[atomIndex] = m_atoms->x[i] - m_systemSize[0]; shouldCopy = true;
    //                }
    
    //                if(shouldCopy) {
    //                    m_atoms->y[atomIndex] = m_atoms->y[i];
    //                    m_atoms->z[atomIndex] = m_atoms->z[i];
    //                }
    //            } else if(dimension==1) {
    //                if(m_atoms->y[i] < m_rShell) {
    //                    m_atoms->y[atomIndex] = m_atoms->y[i] + m_systemSize[1]; shouldCopy = true;
    //                } else if (m_atoms->y[i] > m_systemSize[1] - m_rShell) {
    //                    m_atoms->y[atomIndex] = m_atoms->y[i] - m_systemSize[1]; shouldCopy = true;
    //                }
    
    //                if(shouldCopy) {
    //                    m_atoms->x[atomIndex] = m_atoms->x[i];
    //                    m_atoms->z[atomIndex] = m_atoms->z[i];
    //                }
    //            } else if(dimension==2) {
    //                if(m_atoms->z[i] < m_rShell) {
    //                    m_atoms->z[atomIndex] = m_atoms->z[i] + m_systemSize[2]; shouldCopy = true;
    //                } else if (m_atoms->z[i] > m_systemSize[2] - m_rShell) {
    //                    m_atoms->z[atomIndex] = m_atoms->z[i] - m_systemSize[2]; shouldCopy = true;
    //                }
    
    //                if(shouldCopy) {
    //                    m_atoms->x[atomIndex] = m_atoms->x[i];
    //                    m_atoms->y[atomIndex] = m_atoms->y[i];
    //                }
    //            }
    
    //            if(shouldCopy) {
    //                m_atoms->vx[atomIndex] = m_atoms->vx[i];
    //                m_atoms->vy[atomIndex] = m_atoms->vy[i];
    //                m_atoms->vz[atomIndex] = m_atoms->vz[i];
    
    //                m_atoms->mass[atomIndex] = m_atoms->mass[i];
    //                m_atoms->index[atomIndex] = atomIndex;
    //                m_atoms->numberOfGhostAtoms++;
    //            }
    //        }
    //    }
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    if(m_shouldSample) m_potential->calculateForcesAndEnergyAndPressure(this);
    else m_potential->calculateForces(this);
}

void System::step(float dt) {
    if(!m_initialized) {
        cout << "System not initialized, aborting. Remember to call initialize()." << endl;
        exit(1);
    }
    m_integrator->integrate(this, dt);
    m_steps++;
    m_currentTime += dt;
}
