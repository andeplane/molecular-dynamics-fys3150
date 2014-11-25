#include "neighborlist.h"
#include "celllist.h"
#include "math/vec3.h"
#include "system.h"
#include "cpelapsedtimer.h"
#include "atom.h"
#include <algorithm>
#include <iostream>
#include <cassert>
using namespace std;

NeighborList::NeighborList() :
    m_system(0),
    m_rShellSquared(-1),
    m_numNeighborPairs(0),
    m_neighbors(0)
{

}

void NeighborList::clear() {
    for(unsigned int i=0; i<m_system->atoms().numberOfAtoms; i++) {
        m_neighbors[i].numberOfAtoms = 0;
    }
}

void NeighborList::setup(System *system, float rShell)
{
    m_system = system;
    m_rShellSquared = rShell*rShell;
    m_rCut = system->rCut();
    m_cellList.setup(system, rShell);
    m_neighbors = new MiniAtoms[m_system->atoms().numberOfAtoms];
}

void NeighborList::update()
{
    vec3 systemSize = m_system->systemSize();
    vec3 systemSizeHalf = m_system->systemSize()*0.5;

    // m_system->atoms().sort();
    m_cellList.update();

    CPElapsedTimer::updateNeighborList().start();
    clear();

    Atoms &atoms = m_system->atoms();
    for(int cx=0; cx<m_cellList.numberOfCellsX(); cx++) {
        for(int cy=0; cy<m_cellList.numberOfCellsY(); cy++) {
            for(int cz=0; cz<m_cellList.numberOfCellsZ(); cz++) {
                int cellIndex1 = m_cellList.index(cx, cy, cz);
                const vector<unsigned int> &cell1 = m_cellList.cells()[cellIndex1];

                for(int dx=0; dx<=1; dx++) {
                    for(int dy=(dx==0 ? 0 : -1); dy<=1; dy++) {
                        for(int dz=(dx==0 && dy==0 ? 0 : -1); dz<=1; dz++) {
                            int cellIndex2 = m_cellList.indexPeriodic(cx+dx, cy+dy, cz+dz);
                            const vector<unsigned int> &cell2 = m_cellList.cells()[cellIndex2];

                            const unsigned int cell1Size = cell1.size();
                            for(unsigned int i=0; i<cell1Size; i++) {
                                const unsigned int atom1Index = cell1[i];

                                float x = atoms.x[atom1Index];
                                float y = atoms.y[atom1Index];
                                float z = atoms.z[atom1Index];

                                const unsigned int cell2Size = cell2.size();
#ifdef MD_SIMD
#pragma simd
#endif
                                for(unsigned int j=(dx==0 && dy==0 && dz==0 ? i+1 : 0); j<cell2Size; j++) {
                                    const unsigned int atom2Index = cell2[j];
                                    float dx = x - atoms.x[atom2Index];
                                    float dy = y - atoms.y[atom2Index];
                                    float dz = z - atoms.z[atom2Index];
                                    float ddx = systemSize[0]*( (dx < -systemSizeHalf[0] ) - (dx > systemSizeHalf[0]));
                                    float ddy = systemSize[1]*( (dy < -systemSizeHalf[1] ) - (dy > systemSizeHalf[1]));
                                    float ddz = systemSize[2]*( (dz < -systemSizeHalf[2] ) - (dz > systemSizeHalf[2]));

                                    dx += ddx;
                                    dy += ddy;
                                    dz += ddz;

                                    const float dr2 = dx*dx + dy*dy + dz*dz;

                                    float dr2i = 1.0/dr2;
                                    float dr6i = dr2i*dr2i*dr2i;
                                    MiniAtoms &neighborList1 = m_neighbors[atom1Index];
                                    neighborList1.dx[neighborList1.numberOfAtoms] = dx;
                                    neighborList1.dy[neighborList1.numberOfAtoms] = dy;
                                    neighborList1.dz[neighborList1.numberOfAtoms] = dz;
                                    neighborList1.ddx[neighborList1.numberOfAtoms] = ddx;
                                    neighborList1.ddy[neighborList1.numberOfAtoms] = ddy;
                                    neighborList1.ddz[neighborList1.numberOfAtoms] = ddz;
                                    neighborList1.dr2[neighborList1.numberOfAtoms] = dr2;
                                    neighborList1.dr2i[neighborList1.numberOfAtoms] = dr2i;
                                    neighborList1.dr6i[neighborList1.numberOfAtoms] = dr6i;
                                    neighborList1.originalAtom[neighborList1.numberOfAtoms] = atom2Index;
                                    neighborList1.numberOfAtoms++;

                                    MiniAtoms &neighborList2 = m_neighbors[atom2Index];
                                    neighborList2.dx[neighborList2.numberOfAtoms] = -dx;
                                    neighborList2.dy[neighborList2.numberOfAtoms] = -dy;
                                    neighborList2.dz[neighborList2.numberOfAtoms] = -dz;
                                    neighborList2.ddx[neighborList2.numberOfAtoms] = -ddx;
                                    neighborList2.ddy[neighborList2.numberOfAtoms] = -ddy;
                                    neighborList2.ddz[neighborList2.numberOfAtoms] = -ddz;
                                    neighborList2.dr2[neighborList2.numberOfAtoms] = dr2;
                                    neighborList2.dr2i[neighborList2.numberOfAtoms] = dr2i;
                                    neighborList2.dr6i[neighborList2.numberOfAtoms] = dr6i;
                                    neighborList2.originalAtom[neighborList2.numberOfAtoms] = atom1Index;
                                    neighborList2.numberOfAtoms++;
#ifdef MD_DEBUG
                                    assert(m_neighbors[atom2Index][0] <= MAXNUMNEIGHBORS && "An atom got too many neighbors :/");
#endif
                                }
                                bool sameCell = dx==0 && dy==0 && dz==0;
                                m_numNeighborPairs += cell2Size*(1.0 - 0.5*sameCell);
                            }
                        }}}
            }}}

    CPElapsedTimer::updateNeighborList().stop();
}

void NeighborList::updateCopies()
{
    CPElapsedTimer::updateNeighborCopies().start();
    Atoms &allAtoms = m_system->atoms();
    for(unsigned int i=0; i<allAtoms.numberOfAtoms; i++) {
        MiniAtoms &neighbors = m_neighbors[i];
        for(unsigned int j=0; j<neighbors.numberOfAtoms; j++) {
            const unsigned int originalAtom = neighbors.originalAtom[j];
            float dx = allAtoms.x[i] - allAtoms.x[originalAtom] + neighbors.ddx[j];
            float dy = allAtoms.y[i] - allAtoms.y[originalAtom] + neighbors.ddy[j];
            float dz = allAtoms.z[i] - allAtoms.z[originalAtom] + neighbors.ddz[j];
            float dr2 = dx*dx + dy*dy + dz*dz;
            neighbors.dr2[j] = dr2;
            if(dr2 < m_rShellSquared) {
                float dr2i = 1.0/dr2;
                float dr6i = dr2i*dr2i*dr2i;
                neighbors.dx[j] = dx;
                neighbors.dy[j] = dy;
                neighbors.dz[j] = dz;
                neighbors.dr2i[j] = dr2i;
                neighbors.dr6i[j] = dr6i;
            }
        }
    }
    CPElapsedTimer::updateNeighborCopies().stop();
}

float NeighborList::averageNumNeighbors()
{
    unsigned int numNeighbors = 0;
    for(unsigned int i=0; i<m_system->atoms().numberOfAtoms; i++) {
        numNeighbors += m_neighbors[i].numberOfAtoms;
    }

    return float(numNeighbors)/m_system->atoms().numberOfAtoms;
}
