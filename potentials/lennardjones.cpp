#include "potentials/lennardjones.h"
#include "celllist.h"
#include "neighborlist.h"
#include "cpelapsedtimer.h"
#include <cmath>
#include <iostream>
using namespace std;
LennardJones::LennardJones(MDDataType_t sigma, MDDataType_t epsilon, MDDataType_t cutoffRadius) :
    m_sigma(sigma),
    m_sigma6(pow(sigma, 6.0)),
    m_epsilon(epsilon),
    m_24epsilon(24*epsilon),
    m_rCutSquared(cutoffRadius*cutoffRadius),
    m_timeSinceLastNeighborListUpdate(0),
    m_potentialEnergyAtRcut(0)
{
    MDDataType_t oneOverDrCut2 = 1.0/m_rCutSquared;
    MDDataType_t oneOverDrCut6 = oneOverDrCut2*oneOverDrCut2*oneOverDrCut2;

    m_potentialEnergyAtRcut = 4*m_epsilon*m_sigma6*oneOverDrCut6*(m_sigma6*oneOverDrCut6 - 1);
}

void LennardJones::calculateForcesAllPairs(System *system) {
    CPElapsedTimer::calculateForces().start();
    m_potentialEnergy = 0;
    m_pressureVirial = 0;

    vec3 systemSize = system->systemSize();
    vec3 systemSizeHalf = system->systemSize()*0.5;

    Atoms &atoms = system->atoms();
    const unsigned int numAtoms = atoms.numberOfAtoms;
    for(unsigned int i=0; i<numAtoms; i++) {
        const MDDataType_t x = atoms.x[i];
        const MDDataType_t y = atoms.y[i];
        const MDDataType_t z = atoms.z[i];
        MDDataType_t fix = 0;
        MDDataType_t fiy = 0;
        MDDataType_t fiz = 0;
        MDDataType_t potentialEnergy = 0;
#pragma simd reduction(+: fix, fiy, fiz, potentialEnergy)
        for(unsigned int j=i+1; j<numAtoms; j++) {
            MDDataType_t dx = x - atoms.x[j];
            MDDataType_t dy = y - atoms.y[j];
            MDDataType_t dz = z - atoms.z[j];
            dx += systemSize[0]*( (dx < -systemSizeHalf[0] ) - (dx > systemSizeHalf[0]));
            dy += systemSize[1]*( (dy < -systemSizeHalf[1] ) - (dy > systemSizeHalf[1]));
            dz += systemSize[2]*( (dz < -systemSizeHalf[2] ) - (dz > systemSizeHalf[2]));

            const MDDataType_t dr2 = dx*dx + dy*dy + dz*dz;
            const MDDataType_t oneOverDr2 = 1.0f/dr2;
            const MDDataType_t sigma6OneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2*m_sigma6;
            const MDDataType_t force = -m_24epsilon*sigma6OneOverDr6*(2*sigma6OneOverDr6 - 1.0f)*oneOverDr2*(dr2 < m_rCutSquared);

            fix -= dx*force;
            fiy -= dy*force;
            fiz -= dz*force;

            atoms.fx[j] += dx*force;
            atoms.fy[j] += dy*force;
            atoms.fz[j] += dz*force;

            // m_pressureVirial += force*sqrt(dr2)*dr2;
            potentialEnergy += (4*m_epsilon*sigma6OneOverDr6*(sigma6OneOverDr6 - 1) - m_potentialEnergyAtRcut)*(dr2 < m_rCutSquared);
        }

        atoms.fx[i] += fix;
        atoms.fy[i] += fiy;
        atoms.fz[i] += fiz;
        m_potentialEnergy += potentialEnergy;
        atoms.numberOfComputedForces += numAtoms-i-1;
    }

    CPElapsedTimer::calculateForces().stop();
}

void LennardJones::calculateForces(System *system)
{
#ifdef BENCHMARK
    calculateForcesAllPairs(system);
    return;
#endif
    m_potentialEnergy = 0;
    m_pressureVirial = 0;
    vec3 systemSize = system->systemSize();
    vec3 systemSizeHalf = system->systemSize()*0.5;

    if(!m_timeSinceLastNeighborListUpdate || m_timeSinceLastNeighborListUpdate++ > BUILDNEIGHBORLIST) {
        system->neighborList().update();
        m_timeSinceLastNeighborListUpdate = 1;
    }

    CPElapsedTimer::calculateForces().start();

    Atoms &atoms = system->atoms();
    for(unsigned int i=0; i<atoms.numberOfAtoms; i++) {
        const MDDataType_t x = atoms.x[i];
        const MDDataType_t y = atoms.y[i];
        const MDDataType_t z = atoms.z[i];
        MDDataType_t fix = 0;
        MDDataType_t fiy = 0;
        MDDataType_t fiz = 0;

        const unsigned int *neighbors = system->neighborList().neighborsForAtomWithIndex(i);

        const unsigned int numNeighbors = neighbors[0];
#ifdef MD_SIMD
#pragma simd reduction (+: fix, fiy, fiz)
#endif
        for(unsigned int j=1; j<=numNeighbors; j++) {
            const unsigned int neighborIndex = neighbors[j];
            MDDataType_t dx = x - atoms.x[neighborIndex];
            MDDataType_t dy = y - atoms.y[neighborIndex];
            MDDataType_t dz = z - atoms.z[neighborIndex];

            if(dx < -systemSizeHalf[0]) dx += systemSize[0];
            else if(dx > systemSizeHalf[0]) dx -= systemSize[0];
            if(dy < -systemSizeHalf[1]) dy += systemSize[1];
            else if(dy > systemSizeHalf[1]) dy -= systemSize[1];
            if(dz < -systemSizeHalf[2]) dz += systemSize[2];
            else if(dz > systemSizeHalf[2]) dz -= systemSize[2];

            const MDDataType_t dr2 = dx*dx + dy*dy + dz*dz;
            const MDDataType_t oneOverDr2 = 1.0f/dr2;
            const MDDataType_t sigma6OneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2*m_sigma6;
            const MDDataType_t force = -m_24epsilon*sigma6OneOverDr6*(2*sigma6OneOverDr6 - 1.0f)*oneOverDr2*(dr2 < m_rCutSquared);

            fix -= dx*force;
            fiy -= dy*force;
            fiz -= dz*force;

            atoms.fx[neighborIndex] += dx*force;
            atoms.fy[neighborIndex] += dy*force;
            atoms.fz[neighborIndex] += dz*force;
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
#ifdef BENCHMARK
    calculateForcesAllPairs(system);
    return;
#endif
    m_potentialEnergy = 0;
    m_pressureVirial = 0;
    vec3 systemSize = system->systemSize();
    vec3 systemSizeHalf = system->systemSize()*0.5;

    if(!m_timeSinceLastNeighborListUpdate || m_timeSinceLastNeighborListUpdate++ > BUILDNEIGHBORLIST) {
        system->neighborList().update();
        m_timeSinceLastNeighborListUpdate = 1;
    }

    CPElapsedTimer::calculateForces().start();

    Atoms &atoms = system->atoms();

    for(unsigned i=0; i<atoms.numberOfAtoms; i++) {
        MDDataType_t x = atoms.x[i];
        MDDataType_t y = atoms.y[i];
        MDDataType_t z = atoms.z[i];
        MDDataType_t fix = 0;
        MDDataType_t fiy = 0;
        MDDataType_t fiz = 0;

        unsigned int *neighbors = system->neighborList().neighborsForAtomWithIndex(i);
        MDDataType_t pressureVirial = 0;
        MDDataType_t potentialEnergy = 0;

        const unsigned int numNeighbors = neighbors[0];
#ifdef MD_SIMD
#pragma simd reduction (+: fix, fiy, fiz, pressureVirial, potentialEnergy)
#endif
        for(unsigned int j=1; j<=numNeighbors; j++) {
            unsigned int neighborIndex = neighbors[j];
            MDDataType_t dx = x - atoms.x[neighborIndex];
            MDDataType_t dy = y - atoms.y[neighborIndex];
            MDDataType_t dz = z - atoms.z[neighborIndex];

            if(dx < -systemSizeHalf[0]) dx += systemSize[0];
            else if(dx > systemSizeHalf[0]) dx -= systemSize[0];
            if(dy < -systemSizeHalf[1]) dy += systemSize[1];
            else if(dy > systemSizeHalf[1]) dy -= systemSize[1];
            if(dz < -systemSizeHalf[2]) dz += systemSize[2];
            else if(dz > systemSizeHalf[2]) dz -= systemSize[2];

            const MDDataType_t dr2 = dx*dx + dy*dy + dz*dz;
            const MDDataType_t oneOverDr2 = 1.0f/dr2;
            const MDDataType_t sigma6OneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2*m_sigma6;
            const MDDataType_t force = -m_24epsilon*sigma6OneOverDr6*(2*sigma6OneOverDr6 - 1.0f)*oneOverDr2*(dr2 < m_rCutSquared);

            fix -= dx*force;
            fiy -= dy*force;
            fiz -= dz*force;

            atoms.fx[neighborIndex] += dx*force;
            atoms.fy[neighborIndex] += dy*force;
            atoms.fz[neighborIndex] += dz*force;

            pressureVirial += force*sqrt(dr2)*dr2;
            potentialEnergy += (4*m_epsilon*sigma6OneOverDr6*(sigma6OneOverDr6 - 1.0f) - m_potentialEnergyAtRcut)*(dr2 < m_rCutSquared);
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
