#include "neighborlist.h"
#include "celllist.h"
#include "math/vec3.h"
#include "system.h"
#include "cpelapsedtimer.h"
#include "atom.h"
#include <algorithm>
#include <iostream>
using namespace std;

NeighborList::NeighborList() :
    m_system(0),
    m_rShellSquared(-1)
{
    m_neighbors = new unsigned int*[MAXNUMATOMS];
    for(int i=0; i<MAXNUMATOMS; i++) {
        m_neighbors[i] = new unsigned int[500];
    }
}

void NeighborList::clear() {
    for(unsigned int i=0; i<m_system->atoms().numberOfAtoms; i++) {
        m_neighbors[i][0] = 0;
    }
}

void NeighborList::setup(System *system, float rShell)
{
    m_system = system;
    m_rShellSquared = rShell*rShell;
    m_cellList.setup(system, rShell);
}

void NeighborList::update()
{
    vec3 systemSize = m_system->systemSize();
    vec3 systemSizeHalf = m_system->systemSize()*0.5;

    m_system->atoms().sort();
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
                                unsigned int atom1Index = cell1[i];

                                float x = atoms.x[atom1Index];
                                float y = atoms.y[atom1Index];
                                float z = atoms.z[atom1Index];

                                const unsigned int cell2Size = cell2.size();
#pragma simd
                                for(unsigned int j=(dx==0 && dy==0 && dz==0 ? i+1 : 0); j<cell2Size; j++) {
                                    unsigned int atom2Index = cell2[j];
                                    float dx = x - atoms.x[atom2Index];
                                    float dy = y - atoms.y[atom2Index];
                                    float dz = z - atoms.z[atom2Index];

                                    dx += systemSize[0]*( (dx < -systemSizeHalf[0] ) - (dx > systemSizeHalf[0]));
                                    dy += systemSize[1]*( (dy < -systemSizeHalf[1] ) - (dy > systemSizeHalf[1]));
                                    dz += systemSize[2]*( (dz < -systemSizeHalf[2] ) - (dz > systemSizeHalf[2]));

                                    const float dr2 = dx*dx + dy*dy + dz*dz;

                                    bool shouldNotAdd = dr2 > m_rShellSquared;
                                    m_neighbors[atom2Index][ ++m_neighbors[atom2Index][0] ] = atom1Index;
                                    m_neighbors[atom2Index][0] -= shouldNotAdd; // Decrease neighborcounter if we shouldn't add this
                                }
                            }
                        }}}
            }}}

    CPElapsedTimer::updateNeighborList().stop();
}
