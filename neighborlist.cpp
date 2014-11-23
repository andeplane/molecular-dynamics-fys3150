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
    m_rShellSquared(-1)
{
    m_neighbors = new unsigned int*[MAXNUMATOMS];
    for(int i=0; i<MAXNUMATOMS; i++) {
        m_neighbors[i] = new unsigned int[1000];
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
                    if(cx==m_cellList.numberOfCellsX()-1 && dx==1) continue;
                    for(int dy=(dx==0 ? 0 : -1); dy<=1; dy++) {
                        if( (cy==0 && dy==-1) || (cy==m_cellList.numberOfCellsY()-1 && dy==1) ) continue;
                        for(int dz=(dx==0 && dy==0 ? 0 : -1); dz<=1; dz++) {
                            if( (cz==0 && dz==-1) || (cz==m_cellList.numberOfCellsZ()-1 && dz==1) ) continue;

                            int cellIndex2 = m_cellList.index(cx+dx, cy+dy, cz+dz);
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

                                    const float dr2 = dx*dx + dy*dy + dz*dz;

                                    bool shouldNotAdd = dr2 > m_rShellSquared;
                                    m_neighbors[atom2Index][ ++m_neighbors[atom2Index][0] ] = atom1Index;
                                    m_neighbors[atom2Index][0] -= shouldNotAdd; // Decrease neighborcounter if we shouldn't add this
                                    assert(m_neighbors[atom2Index][0] <= 500 && "Too many neighbors for one atom, sorry.");
                                }
                            }
                        }}}
            }}}

    CPElapsedTimer::updateNeighborList().stop();
}
