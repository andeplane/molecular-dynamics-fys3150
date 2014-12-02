#include "celllist.h"
#include "system.h"
#include <iostream>
#include "cpelapsedtimer.h"
#include <cassert>

using namespace std;

unsigned int CellList::indexf(MDDataType_t x, MDDataType_t y, MDDataType_t z)
{
    unsigned int cx = x*m_oneOverLengthX;
    unsigned int cy = y*m_oneOverLengthY;
    unsigned int cz = z*m_oneOverLengthZ;

    return index(cx, cy, cz);
}

vector<vector<unsigned int> > &CellList::cells()
{
    return m_cells;
}

CellList::CellList() :
    m_system(0),
    m_numberOfCellsX(-1),
    m_numberOfCellsY(-1),
    m_numberOfCellsZ(-1),
    m_oneOverLengthX(0),
    m_oneOverLengthY(0),
    m_oneOverLengthZ(0)
{

}

void CellList::setup(System *system, MDDataType_t rCut)
{
    m_system = system;

    m_numberOfCellsX = system->systemSize().x() / rCut;
    m_numberOfCellsY = system->systemSize().y() / rCut;
    m_numberOfCellsZ = system->systemSize().z() / rCut;

    m_oneOverLengthX = 1.0 / (system->systemSize().x()/m_numberOfCellsX);
    m_oneOverLengthY = 1.0 / (system->systemSize().y()/m_numberOfCellsY);
    m_oneOverLengthZ = 1.0 / (system->systemSize().z()/m_numberOfCellsZ);

    m_cells.resize(m_numberOfCellsX*m_numberOfCellsY*m_numberOfCellsZ);
    m_neighbors.resize(m_cells.size());

    for(int cx=0; cx<numberOfCellsX(); cx++) {
        for(int cy=0; cy<numberOfCellsY(); cy++) {
            for(int cz=0; cz<numberOfCellsZ(); cz++) {
                int cellIndex1 = index(cx, cy, cz);

                for(int dx=0; dx<=1; dx++) {
                    for(int dy=(dx==0 ? 0 : -1); dy<=1; dy++) {
                        for(int dz=(dx==0 && dy==0 ? 0 : -1); dz<=1; dz++) {
                            int cellIndex2 = indexPeriodic(cx+dx, cy+dy, cz+dz);
                            vector<unsigned int> *cell2 = &m_cells[cellIndex2];
                            m_neighbors[cellIndex1].push_back(cell2);
                        }
                    }
                }
                assert(m_neighbors[cellIndex1].size() == 14 && "Balle");
            }
        }
    }
}

void CellList::clear()
{
    for(unsigned int i=0; i<m_cells.size(); i++) {
        m_cells[i].clear();
    }
}

void CellList::update()
{
    CPElapsedTimer::updateCellList().start();
    clear();
    Atoms &atoms = m_system->atoms();
    for(unsigned int i=0; i<atoms.numberOfAtoms; i++) {
        unsigned int cellIndex = indexf(atoms.x[i], atoms.y[i], atoms.z[i]);
        atoms.cellIndex[i] = cellIndex;
#ifdef MD_DEBUG
        if(!(cellIndex < m_cells.size() && cellIndex >= 0)) {
            cout << "Atom " << i << " with position [" << atoms.x[i] << ", " << atoms.y[i] << ", " << atoms.z[i] << "] has cell index " << cellIndex << " which has bounds [0," << m_cells.size() << "]." << endl;
            cout << "System size is: " << m_system->systemSize() << endl;
            printf("Super detailed position: [%.16f, %.16f, %.16f]\n", atoms.x[i], atoms.y[i], atoms.z[i]);
            printf("Super system size: [%.16f, %.16f, %.16f]\n", m_system->systemSize().x(), m_system->systemSize().y(), m_system->systemSize().z());
            cout << "Smaller than test would say: [" << (atoms.x[i] < m_system->systemSize().x()) << ", " << (atoms.y[i] < m_system->systemSize().y()) << ", " << (atoms.z[i] < m_system->systemSize().z()) << "]" << endl;
            cout << "Larger than or equal test would say: [" << (atoms.x[i] >= m_system->systemSize().x()) << ", " << (atoms.y[i] >= m_system->systemSize().y()) << ", " << (atoms.z[i] >= m_system->systemSize().z()) << "]" << endl;
        }
        assert(cellIndex < m_cells.size() && cellIndex >= 0 && "Cell index out of bounds");
#endif
        m_cells[cellIndex].push_back(i);
    }

    CPElapsedTimer::updateCellList().stop();
}
