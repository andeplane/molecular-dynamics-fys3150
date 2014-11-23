#include "celllist.h"
#include "atom.h"
#include "system.h"
#include <iostream>
#include "cpelapsedtimer.h"
#include <cassert>

using namespace std;

unsigned int CellList::indexf(float x, float y, float z)
{
    unsigned int cx = x/m_system->systemSize()[0]*(m_numberOfCellsX-2) + 1;
    unsigned int cy = y/m_system->systemSize()[1]*(m_numberOfCellsY-2) + 1;
    unsigned int cz = z/m_system->systemSize()[2]*(m_numberOfCellsZ-2) + 1;

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
    m_numberOfCellsZ(-1)
{

}

void CellList::setup(System *system, float minimumCellSize)
{
    m_system = system;

    m_numberOfCellsX = system->systemSize().x() / minimumCellSize + 2;
    m_numberOfCellsY = system->systemSize().y() / minimumCellSize + 2;
    m_numberOfCellsZ = system->systemSize().z() / minimumCellSize + 2;

    m_cells.resize(m_numberOfCellsX*m_numberOfCellsY*m_numberOfCellsZ);
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
    // cout << "System size: " << m_system->systemSize() << " and rShell = " << m_system->m_rShell << endl;
    Atoms &atoms = m_system->atoms();
    for(unsigned int i=0; i<atoms.numberOfAtomsIncludingGhosts(); i++) {
        unsigned int cellIndex = indexf(atoms.x[i], atoms.y[i], atoms.z[i]);
        atoms.cellIndex[i] = cellIndex;
        // cout << "Atom " << i << " at " << m_system->atoms().x[i] << ", " << m_system->atoms().y[i] << ", " << m_system->atoms().z[i] << " has cell index " << cellIndex << " while we have " << m_cells.size() << " cells." << endl;
        assert(cellIndex < m_cells.size() && cellIndex >= 0 && "Cell index out of bounds");
        m_cells[cellIndex].push_back(i);
    }

    CPElapsedTimer::updateCellList().stop();
}
