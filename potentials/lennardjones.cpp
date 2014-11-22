#include "potentials/lennardjones.h"
#include "celllist.h"
#include "neighborlist.h"
#include "cpelapsedtimer.h"
#include <cmath>
#include <iostream>
using namespace std;
LennardJones::LennardJones(float sigma, float epsilon, float cutoffRadius) :
    m_sigma(sigma),
    m_sigma6(pow(sigma, 6.0)),
    m_epsilon(epsilon),
    m_24epsilon(24*epsilon),
    m_rCutSquared(cutoffRadius*cutoffRadius),
    m_timeSinceLastNeighborListUpdate(0),
    m_potentialEnergyAtRcut(0)
{
    float oneOverDrCut2 = 1.0/m_rCutSquared;
    float oneOverDrCut6 = oneOverDrCut2*oneOverDrCut2*oneOverDrCut2;

    m_potentialEnergyAtRcut = 4*m_epsilon*m_sigma6*oneOverDrCut6*(m_sigma6*oneOverDrCut6 - 1);
}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0;
    m_pressureVirial = 0;
    vec3 systemSize = system->systemSize();
    vec3 systemSizeHalf = system->systemSize()*0.5;

    if(!m_timeSinceLastNeighborListUpdate || m_timeSinceLastNeighborListUpdate++ > 20) {
        system->neighborList().update();
        m_timeSinceLastNeighborListUpdate = 1;
    }

    CPElapsedTimer::calculateForces().start();
    const unsigned int numAtoms = system->atoms().size();
    for(unsigned int i=0; i<numAtoms; i++) {
        Atom *atom1 = system->atoms()[i];

        const vector<Atom*> &neighbors = system->neighborList().neighborsForAtomWithIndex(atom1->index());
        const unsigned int numNeighbors = neighbors.size();
        for(unsigned int j=0; j<numNeighbors; j++) {
            Atom *atom2 = neighbors[j];

            vec3 deltaRVector = atom1->position;
            deltaRVector.addAndMultiply(atom2->position, -1);

            // Minimum image convention
            if(deltaRVector[0] > systemSizeHalf[0]) deltaRVector[0] -= systemSize[0];
            else if(deltaRVector[0] < -systemSizeHalf[0]) deltaRVector[0] += systemSize[0];
            if(deltaRVector[1] > systemSizeHalf[1]) deltaRVector[1] -= systemSize[1];
            else if(deltaRVector[1] < -systemSizeHalf[1]) deltaRVector[1] += systemSize[1];
            if(deltaRVector[2] > systemSizeHalf[2]) deltaRVector[2] -= systemSize[2];
            else if(deltaRVector[2] < -systemSizeHalf[2]) deltaRVector[2] += systemSize[2];

            const float dr2 = deltaRVector.lengthSquared();
            const float oneOverDr2 = 1.0f/dr2;
            const float oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;

            const float force = -m_24epsilon*m_sigma6*oneOverDr6*(2*m_sigma6*oneOverDr6 - 1)*oneOverDr2*(dr2 < m_rCutSquared);

            atom1->force.addAndMultiply(deltaRVector, -force);
            atom2->force.addAndMultiply(deltaRVector, force);

            if(m_shouldComputeEnergyAndPressureVirial) {
                m_pressureVirial += force*sqrt(dr2)*dr2;
                m_potentialEnergy += (4*m_epsilon*m_sigma6*oneOverDr6*(m_sigma6*oneOverDr6 - 1) - m_potentialEnergyAtRcut)*(dr2 < m_rCutSquared);
            }

        }
    }

    CPElapsedTimer::calculateForces().stop();
}
