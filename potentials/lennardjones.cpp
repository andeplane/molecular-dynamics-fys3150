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
    m_rCutSquared(cutoffRadius*cutoffRadius),
    m_timeSinceLastNeighborListUpdate(0),
    m_potentialEnergyAtRcut(0)
{
    float oneOverDrCut2 = 1.0/m_rCutSquared;
    float oneOverDrCut6 = oneOverDrCut2*oneOverDrCut2*oneOverDrCut2;

    m_potentialEnergyAtRcut = 4*m_epsilon*m_sigma6*oneOverDrCut6*(m_sigma6*oneOverDrCut6 - 1);
}

void LennardJones::calculateForces(System *system) {
    int size = system->atoms().size();

    for(int i=0; i<size; i++) {
        vec3 systemSize = system->systemSize();

#pragma ivdep
        for(int j=0; j<size; j++) {
            float dx = system->m_positionsAndForces[6*i+0] - system->m_positionsAndForces[6*j+0];
            float dy = system->m_positionsAndForces[6*i+1] - system->m_positionsAndForces[6*j+1];
            float dz = system->m_positionsAndForces[6*i+2] - system->m_positionsAndForces[6*j+2];

            dx += systemSize[0]*((dx < -0.5*systemSize[0]) - (dx > 0.5*systemSize[0]));
            dy += systemSize[1]*((dy < -0.5*systemSize[1]) - (dy > 0.5*systemSize[1]));
            dz += systemSize[2]*((dz < -0.5*systemSize[2]) - (dz > 0.5*systemSize[2]));

            float dr2 = dx*dx + dy*dy + dz*dz + 1e-9;
            float oneOverDr2 = 1.0/dr2*(i!=j);
            float oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;

            float force = -24*m_epsilon*m_sigma6*oneOverDr6*(2*m_sigma6*oneOverDr6 - 1)*oneOverDr2 * (dr2 < m_rCutSquared);

            system->m_positionsAndForces[6*j+3] += force*dx;
            system->m_positionsAndForces[6*j+4] += force*dy;
            system->m_positionsAndForces[6*j+5] += force*dz;
        }
    }
}

//void LennardJones::calculateForces2(System *system)
//{
//    m_potentialEnergy = 0;
//    m_pressureVirial = 0;
//    vec3 systemSize = system->systemSize();

//    if(!m_timeSinceLastNeighborListUpdate || m_timeSinceLastNeighborListUpdate++ > 20) {
//        system->neighborList().update();
//        m_timeSinceLastNeighborListUpdate = 1;
//    }

//    CPElapsedTimer::calculateForces().start();
//    for(unsigned int i=0; i<system->atoms().size(); i++) {
//        Atom *atom1 = system->atoms()[i];

//        vector<float *> &neighbors = system->neighborList().neighborsForAtomWithIndex(atom1->index());

//        int size = neighbors.size();
//        #pragma ivdep
//        for(unsigned int k=0; k<size; k++) {
//            float *pos = neighbors[k];

//            float dx = system->m_positionsAndForces[6*i+0] - pos[0];
//            float dy = system->m_positionsAndForces[6*i+1] - pos[1];
//            float dz = system->m_positionsAndForces[6*i+2] - pos[2];

//            dx += systemSize[0]*((dx < -0.5*systemSize[0]) - (dx > 0.5*systemSize[0]));
//            dy += systemSize[1]*((dy < -0.5*systemSize[1]) - (dy > 0.5*systemSize[1]));
//            dz += systemSize[2]*((dz < -0.5*systemSize[2]) - (dz > 0.5*systemSize[2]));

//            float dr2 = dx*dx + dy*dy + dz*dz;
//            float oneOverDr2 = 1.0/dr2;
//            float oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;

//            float force = -24*m_epsilon*m_sigma6*oneOverDr6*(2*m_sigma6*oneOverDr6 - 1)*oneOverDr2 * (dr2 < m_rCutSquared);

//            pos[3] += force*dx;
//            pos[4] += force*dy;
//            pos[5] += force*dz;

//            // m_pressureVirial += force*sqrt(dr2)*dr2;
//            // m_potentialEnergy += (4*m_epsilon*m_sigma6*oneOverDr6*(m_sigma6*oneOverDr6 - 1) - m_potentialEnergyAtRcut)*(dr2 < m_rCutSquared);
//        }
//    }

//    CPElapsedTimer::calculateForces().stop();
//}
