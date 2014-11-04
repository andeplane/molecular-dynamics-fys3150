#include "neighborlist.h"
#include "celllist.h"
#include "math/vec3.h"
#include "atom.h"
#include "system.h"
#include <iostream>
using namespace std;

NeighborList::NeighborList() :
    m_system(0),
    m_rShellSquared(-1)
{

}

void NeighborList::clear() {
    for(int i=0; i<m_neighbors.size(); i++) {
        m_neighbors[i].clear();
    }
}

void NeighborList::setup(System *system, float rShell)
{
    m_system = system;
    m_rShellSquared = rShell*rShell;
    m_cellList.setup(system, rShell);
    m_neighbors.resize(system->atoms().size());
    update();
}

void NeighborList::update()
{
    vec3 systemSize = m_system->systemSize();

    m_cellList.update();
    clear();

    for(int cx=0; cx<m_cellList.numberOfCellsX(); cx++) {
    for(int cy=0; cy<m_cellList.numberOfCellsY(); cy++) {
    for(int cz=0; cz<m_cellList.numberOfCellsZ(); cz++) {
        int cellIndex1 = m_cellList.index(cx, cy, cz);
        vector<Atom*> &cell1 = m_cellList.cells().at(cellIndex1);

        for(int dx=-1; dx<=1; dx++) {
        for(int dy=-1; dy<=1; dy++) {
        for(int dz=-1; dz<=1; dz++) {
            int cellIndex2 = m_cellList.indexPeriodic(cx+dx, cy+dy, cz+dz);
            vector<Atom*> &cell2 = m_cellList.cells().at(cellIndex2);

            for(int i=0; i<cell1.size(); i++) {
                Atom *atom1 = cell1[i];
                for(int j=0; j<cell2.size(); j++) {
                    Atom *atom2 = cell2[j];
                    if(atom1->index() <= atom2->index()) continue;

                    vec3 deltaRVector = atom1->position - atom2->position;

                    for(int a=0; a<3; a++) {
                        if(deltaRVector[a] > 0.5*systemSize[a]) deltaRVector[a] -= systemSize[a];
                        else if(deltaRVector[a] < -0.5*systemSize[a]) deltaRVector[a] += systemSize[a];
                    }

                    float dr2 = deltaRVector.lengthSquared();
                    if(dr2 > m_rShellSquared) continue;
                    m_neighbors[atom1->index()].push_back(atom2);
                }
            }
        }}}
    }}}

//    int count = 0;
//    for(int i=0; i<m_neighbors.size(); i++) {
//        count += m_neighbors[i].size();
//    }

//    count /= m_neighbors.size();
//    cout << "We had an average of " << count << " neighbors per atom." << endl;
}
