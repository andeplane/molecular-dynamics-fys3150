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
    m_numNeighborPairs(0)
{
    m_neighbors = new unsigned int*[MAXNUMATOMS];
    for(int i=0; i<MAXNUMATOMS; i++) {
        m_neighbors[i] = new unsigned int[MAXNUMNEIGHBORS];
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

    // m_system->atoms().sort();
    m_cellList.update();

    CPElapsedTimer::updateNeighborList().start();
    clear();

    Atoms &atoms = m_system->atoms();
    const unsigned int cellSize = m_cellList.cells().size();
    for(unsigned int cellIndex1=0; cellIndex1<cellSize; cellIndex1++) {
        const vector<unsigned int> &cell1 = m_cellList.cells()[cellIndex1];
        vector<vector<unsigned int> *> &neighbors = m_cellList.getNeighbors(cellIndex1);
        for(auto cell2Pointer : neighbors) {
            vector<unsigned int> &cell2 = *cell2Pointer;
            const unsigned int cell1Size = cell1.size();
            for(unsigned int i=0; i<cell1Size; i++) {
                const unsigned int atom1Index = cell1[i];

                float x = atoms.x[atom1Index];
                float y = atoms.y[atom1Index];
                float z = atoms.z[atom1Index];
                const bool sameCell = &cell1==&cell2;
                const unsigned int cell2Size = cell2.size();
#ifdef MD_SIMD
#pragma simd
#endif
                for(unsigned int j=(sameCell ? i+1 : 0); j<cell2Size; j++) {
                    const unsigned int atom2Index = cell2[j];
                    float dx = x - atoms.x[atom2Index];
                    float dy = y - atoms.y[atom2Index];
                    float dz = z - atoms.z[atom2Index];

                    dx += systemSize[0]*( (dx < -systemSizeHalf[0] ) - (dx > systemSizeHalf[0]));
                    dy += systemSize[1]*( (dy < -systemSizeHalf[1] ) - (dy > systemSizeHalf[1]));
                    dz += systemSize[2]*( (dz < -systemSizeHalf[2] ) - (dz > systemSizeHalf[2]));

                    const float dr2 = dx*dx + dy*dy + dz*dz;

                    const bool shouldNotAdd = dr2 > m_rShellSquared;
                    m_neighbors[atom2Index][ ++m_neighbors[atom2Index][0] ] = atom1Index;
                    m_neighbors[atom2Index][0] -= shouldNotAdd; // Decrease neighborcounter if we shouldn't add this
#ifdef MD_DEBUG
                    assert(m_neighbors[atom2Index][0] <= MAXNUMNEIGHBORS && "An atom got too many neighbors :/");
#endif
                }
                m_numNeighborPairs += cell2Size*(1.0 - 0.5*sameCell);
            }
        }
    }
    CPElapsedTimer::updateNeighborList().stop();
}

float NeighborList::averageNumNeighbors()
{
    unsigned int numNeighbors = 0;
    for(unsigned int i=0; i<m_system->atoms().numberOfAtoms; i++) {
        numNeighbors += m_neighbors[i][0];
    }

    return float(numNeighbors)/m_system->atoms().numberOfAtoms;
}
