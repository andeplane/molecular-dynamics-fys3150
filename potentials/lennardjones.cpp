#include "potentials/lennardjones.h"
#include "celllist.h"
#include "neighborlist.h"
#include <cmath>
#include <iostream>
using namespace std;
LennardJones::LennardJones(float sigma, float epsilon, float cutoffRadius) :
    m_sigma(sigma),
    m_epsilon(epsilon),
    m_rCutSquared(cutoffRadius*cutoffRadius),
    m_timeSinceLastNeighborListUpdate(0)
{

}

#define CELLLISTS
// #define NEIGHBORLISTS

inline void LennardJones::calculateForcesBetweenAtoms(Atom *atom1, Atom *atom2, const vec3 &systemSize) {
    vec3 deltaRVector = atom1->position - atom2->position;

    // Minimum image convention
    for(int a=0; a<3; a++) {
        if(deltaRVector[a] > 0.5*systemSize[a]) deltaRVector[a] -= systemSize[a];
        else if(deltaRVector[a] < -0.5*systemSize[a]) deltaRVector[a] += systemSize[a];
    }

    float dr2 = deltaRVector.lengthSquared();
    if(dr2 > m_rCutSquared) return;

    float oneOverDr2 = 1.0/dr2;
    float oneOverDr6 = oneOverDr2*oneOverDr2*oneOverDr2;

    float sigmaSixth = pow(m_sigma, 6.0);
    float force = -24*m_epsilon*sigmaSixth*oneOverDr6*(2*sigmaSixth*oneOverDr6 - 1)*oneOverDr2;

    atom1->force.addAndMultiply(deltaRVector, -force);
    atom2->force.addAndMultiply(deltaRVector, force);

    float potentialEnergy = 4*m_epsilon*sigmaSixth*oneOverDr6*(sigmaSixth*oneOverDr6 - 1);
    addPotentialEnergy(potentialEnergy);
}


#ifdef CELLLISTS
void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    CellList &cellList = system->cellList();
    cellList.update();
    system->resetForcesOnAllAtoms();

    for(int cx=0; cx<cellList.numberOfCellsX(); cx++) {
        for(int cy=0; cy<cellList.numberOfCellsY(); cy++) {
            for(int cz=0; cz<cellList.numberOfCellsZ(); cz++) {
                int cellIndex1 = cellList.index(cx, cy, cz);
                vector<Atom*> &cell1 = cellList.cells().at(cellIndex1);
                for(int dx=-1; dx<=1; dx++) {
                    for(int dy=-1; dy<=1; dy++) {
                        for(int dz=-1; dz<=1; dz++) {
                            int cellIndex2 = cellList.indexPeriodic(cx+dx, cy+dy, cz+dz);
                            vector<Atom*> &cell2 = cellList.cells().at(cellIndex2);
                            vec3 systemSize = system->systemSize();

                            for(int i=0; i<cell1.size(); i++) {
                                Atom *atom1 = cell1[i];
                                for(int j=0; j<cell2.size(); j++) {
                                    Atom *atom2 = cell2[j];
                                    if(atom1->index() <= atom2->index()) continue; // Newton's 3rd law
                                    calculateForcesBetweenAtoms(atom1,atom2, systemSize);

                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
#else
#ifdef NEIGHBORLISTS

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    vec3 systemSize = system->systemSize();
    int count = 0;
    if(m_timeSinceLastNeighborListUpdate++ > 20) {
        system->neighborList().update();
        m_timeSinceLastNeighborListUpdate = 0;
    }

    for(int i=0; i<system->atoms().size(); i++) {
        Atom *atom1 = system->atoms()[i];
        vector<Atom*> &neighbors = system->neighborList().neighborsForAtomWithIndex(atom1->index());
        for(int j=0; j<neighbors.size(); j++) {
            Atom *atom2 = neighbors[j];
            calculateForcesBetweenAtoms(atom1, atom2, systemSize);
        }
    }
}

#else
void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    vec3 systemSize = system->systemSize();
    int count = 0;
    for(int i=0; i<system->atoms().size(); i++) {
        for(int j=i+1; j<system->atoms().size(); j++) {
            Atom *atom1 = system->atoms()[i];
            Atom *atom2 = system->atoms()[j];
            calculateForcesBetweenAtoms(atom1, atom2, systemSize);
            count++;
        }
    }
}
#endif
#endif
