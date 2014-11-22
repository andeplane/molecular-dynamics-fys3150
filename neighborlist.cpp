#include "neighborlist.h"
#include "celllist.h"
#include "math/vec3.h"
#include "atom.h"
#include "system.h"
#include "cpelapsedtimer.h"

#include <iostream>
using namespace std;

NeighborList::NeighborList() :
    m_system(0),
    m_rShellSquared(-1)
{

}

void NeighborList::clear() {
    for(unsigned int i=0; i<m_neighbors.size(); i++) {
        m_neighbors[i].clear();
    }
}

void NeighborList::setup(System *system, float rShell)
{
    m_system = system;
    m_rShellSquared = rShell*rShell;
    m_cellList.setup(system, rShell);
    m_neighbors.resize(system->atoms().size());
}

void NeighborList::update()
{
    vec3 systemSize = m_system->systemSize();

    m_cellList.update();

    CPElapsedTimer::updateNeighborList().start();
    clear();

//    for(int i=0; i<m_system->atoms().size(); i++) {
//        Atom *atom = m_system->atoms()[i];
//        atom->resetNeighbors();
//    }

    for(int cx=0; cx<m_cellList.numberOfCellsX(); cx++) {
    for(int cy=0; cy<m_cellList.numberOfCellsY(); cy++) {
    for(int cz=0; cz<m_cellList.numberOfCellsZ(); cz++) {
        int cellIndex1 = m_cellList.index(cx, cy, cz);
        const vector<Atom*> &cell1 = m_cellList.cells()[cellIndex1];

        for(int dx=-1; dx<=1; dx++) {
        for(int dy=-1; dy<=1; dy++) {
        for(int dz=-1; dz<=1; dz++) {
            int cellIndex2 = m_cellList.indexPeriodic(cx+dx, cy+dy, cz+dz);
            const vector<Atom*> &cell2 = m_cellList.cells()[cellIndex2];
            unsigned int cell1Size = cell1.size();
            for(unsigned int i=0; i<cell1Size; i++) {
                Atom *atom1 = cell1[i];
                unsigned int cell2Size = cell2.size();

                for(unsigned int j=0; j<cell2Size; j++) {
                    Atom *atom2 = cell2[j];
                    if(atom1->index() <= atom2->index()) continue;

                    vec3 deltaRVector = atom1->position;
                    deltaRVector.addAndMultiply(atom2->position, -1);

                    if(deltaRVector[0] > 0.5*systemSize[0]) deltaRVector[0] -= systemSize[0];
                    else if(deltaRVector[0] < -0.5*systemSize[0]) deltaRVector[0] += systemSize[0];
                    if(deltaRVector[1] > 0.5*systemSize[1]) deltaRVector[1] -= systemSize[1];
                    else if(deltaRVector[1] < -0.5*systemSize[1]) deltaRVector[1] += systemSize[1];
                    if(deltaRVector[2] > 0.5*systemSize[2]) deltaRVector[2] -= systemSize[2];
                    else if(deltaRVector[2] < -0.5*systemSize[2]) deltaRVector[2] += systemSize[2];

                    float dr2 = deltaRVector.lengthSquared();
                    if(dr2 > m_rShellSquared) continue;
                    // atom1->addNeighbor(atom2);
                    m_neighbors[atom1->index()].push_back(atom2);
                }
            }
        }}}
    }}}

    CPElapsedTimer::updateNeighborList().stop();
}
