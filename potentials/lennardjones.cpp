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

void LennardJones::calculateForcesAllPairs(System *system) {
    m_potentialEnergy = 0;
    m_pressureVirial = 0;

    vec3 systemSize = system->systemSize();
    vec3 systemSizeHalf = system->systemSize()*0.5;

    Atoms &atoms = system->atoms();
    const unsigned int numAtoms = atoms.numberOfAtoms;
    for(unsigned int i=0; i<numAtoms; i++) {
        const float x = atoms.x[i];
        const float y = atoms.y[i];
        const float z = atoms.z[i];
        float fix = 0;
        float fiy = 0;
        float fiz = 0;
#pragma simd reduction(+: fix, fiy, fiz)
        for(unsigned int j=i+1; j<numAtoms; j++) {
            float dx = x - atoms.x[j];
            float dy = y - atoms.y[j];
            float dz = z - atoms.z[j];
            dx += systemSize[0]*( (dx < -systemSizeHalf[0] ) - (dx > systemSizeHalf[0]));
            dy += systemSize[1]*( (dy < -systemSizeHalf[1] ) - (dy > systemSizeHalf[1]));
            dz += systemSize[2]*( (dz < -systemSizeHalf[2] ) - (dz > systemSizeHalf[2]));

            const float dr2 = dx*dx + dy*dy + dz*dz;
            const float oneOverDr2 = 1.0f/dr2;
            const float sigma6OneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2*m_sigma6;
            const float force = -m_24epsilon*sigma6OneOverDr6*(2*sigma6OneOverDr6 - 1.0f)*oneOverDr2*(dr2 < m_rCutSquared);

            fix -= dx*force;
            fiy -= dy*force;
            fiz -= dz*force;

            atoms.fx[j] += dx*force;
            atoms.fy[j] += dy*force;
            atoms.fz[j] += dz*force;

            // m_pressureVirial += force*sqrt(dr2)*dr2;
            // m_potentialEnergy += (4*m_epsilon*sigma6OneOverDr6*(sigma6OneOverDr6 - 1) - m_potentialEnergyAtRcut)*(dr2 < m_rCutSquared);
        }

        atoms.fx[i] += fix;
        atoms.fy[i] += fiy;
        atoms.fz[i] += fiz;
        atoms.numberOfComputedForces += numAtoms-i-1;
    }
}

void LennardJones::calculateForces(System *system)
{
//    calculateForcesAllPairs(system);
//    return;
    m_potentialEnergy = 0;
    m_pressureVirial = 0;

    if(!m_timeSinceLastNeighborListUpdate || m_timeSinceLastNeighborListUpdate++ > 20) {
        system->neighborList().update();
        m_timeSinceLastNeighborListUpdate = 1;
    } else {
        system->neighborList().updateCopies();
    }
    return;

    CPElapsedTimer::calculateForces().start();

    Atoms &atoms = system->atoms();
    for(unsigned int i=0; i<atoms.numberOfAtoms; i++) {
        const float x = atoms.x[i];
        const float y = atoms.y[i];
        const float z = atoms.z[i];
        float fix = 0;
        float fiy = 0;
        float fiz = 0;

        const MiniAtoms &neighbors = system->neighborList().neighborsForAtomWithIndex(i);

        const unsigned int numNeighbors = neighbors.numberOfAtoms;
#ifdef MD_SIMD
#pragma simd reduction (+: fix, fiy, fiz)
#endif
        for(unsigned int j=0; j<numNeighbors; j++) {

            float sigma6OneOverDr6 = neighbors.dr6i[j]*m_sigma6;
            const float force = m_24epsilon*sigma6OneOverDr6*(2*sigma6OneOverDr6 - 1.0f)*neighbors.dr2i[j]*(neighbors.dr2[j] < m_rCutSquared);

            fix += neighbors.dx[j]*force;
            fiy += neighbors.dy[j]*force;
            fiz += neighbors.dz[j]*force;
        }

        atoms.numberOfComputedForces += numNeighbors;
        atoms.fx[i] += fix;
        atoms.fy[i] += fiy;
        atoms.fz[i] += fiz;
    }

    CPElapsedTimer::calculateForces().stop();
}

void LennardJones::calculateForcesAndEnergyAndPressure(System *system)
{
    m_potentialEnergy = 0;
    m_pressureVirial = 0;

    if(!m_timeSinceLastNeighborListUpdate || m_timeSinceLastNeighborListUpdate++ > 20) {
        system->neighborList().update();
        m_timeSinceLastNeighborListUpdate = 1;
    } else {
        system->neighborList().updateCopies();
    }

    return;

    CPElapsedTimer::calculateForces().start();

    Atoms &atoms = system->atoms();

    for(unsigned i=0; i<atoms.numberOfAtoms; i++) {
        float fix = 0;
        float fiy = 0;
        float fiz = 0;

        const MiniAtoms &neighbors = system->neighborList().neighborsForAtomWithIndex(i);
        const unsigned int numNeighbors = neighbors.numberOfAtoms;
        float pressureVirial = 0;
        float potentialEnergy = 0;

#ifdef MD_SIMD
// #pragma simd reduction (+: fix, fiy, fiz, pressureVirial, potentialEnergy)
#endif
        for(unsigned int j=0; j<numNeighbors; j++) {
            const float sigma6OneOverDr6 = neighbors.dr6i[j]*m_sigma6;
            const float force = m_24epsilon*sigma6OneOverDr6*(2*sigma6OneOverDr6 - 1.0f)*neighbors.dr2i[j]*(neighbors.dr2[j] < m_rCutSquared);

            fix += neighbors.dx[j]*force;
            fiy += neighbors.dy[j]*force;
            fiz += neighbors.dz[j]*force;

            pressureVirial += 0.5*force*sqrt(neighbors.dr2[j])*neighbors.dr2[j];
            potentialEnergy += 0.5*(4*m_epsilon*sigma6OneOverDr6*(sigma6OneOverDr6 - 1.0f) - m_potentialEnergyAtRcut)*(neighbors.dr2[j] < m_rCutSquared);
        }

        atoms.numberOfComputedForces += numNeighbors;
        atoms.fx[i] += fix;
        atoms.fy[i] += fiy;
        atoms.fz[i] += fiz;
        m_pressureVirial += pressureVirial;
        m_potentialEnergy += potentialEnergy;
    }

    CPElapsedTimer::calculateForces().stop();
}
