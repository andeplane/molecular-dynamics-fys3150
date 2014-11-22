#include "celllist.h"
#include "atom.h"
#include "system.h"
#include <iostream>
#include "cpelapsedtimer.h"

using namespace std;


int CellList::index(const vec3 &position)
{
    int cx = position[0]/m_system->systemSize()[0]*m_numberOfCellsX;
    int cy = position[1]/m_system->systemSize()[1]*m_numberOfCellsY;
    int cz = position[2]/m_system->systemSize()[2]*m_numberOfCellsZ;

    return index(cx, cy, cz);
}

void CellList::sort()
{

}

vector<vector<Atom *> > &CellList::cells()
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

void CellList::setup(System *system, float rCut)
{
    m_system = system;

    m_numberOfCellsX = system->systemSize().x() / rCut;
    m_numberOfCellsY = system->systemSize().y() / rCut;
    m_numberOfCellsZ = system->systemSize().z() / rCut;

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
    for(unsigned int i=0; i<m_system->atoms().size(); i++) {
        Atom &atom = m_system->atoms()[i];
        int cellIndex = index(atom.position);
        atom.setCellIndex(cellIndex);
        m_cells[cellIndex].push_back(&atom);
    }

    CPElapsedTimer::updateCellList().stop();
}
