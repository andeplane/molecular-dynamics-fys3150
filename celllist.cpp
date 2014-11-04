#include "celllist.h"
#include "atom.h"
#include "system.h"
#include <iostream>
#include "cpelapsedtimer.h"

using namespace std;

int CellList::index(int cx, int cy, int cz)
{
    return (cx+1)*m_numberOfCellsY*m_numberOfCellsZ + (cy+1)*m_numberOfCellsZ + (cz+1);
    // return cx*m_numberOfCellsY*m_numberOfCellsZ + cy*m_numberOfCellsZ + cz;
}

int CellList::indexPeriodic(int cx, int cy, int cz)
{
    cout << "Should not be used." << endl;
    exit(1);
    // return ( (cx+m_numberOfCellsX) % m_numberOfCellsX)*m_numberOfCellsY*m_numberOfCellsZ + ( (cy+m_numberOfCellsY) % m_numberOfCellsY)*m_numberOfCellsZ + ( (cz+m_numberOfCellsZ) % m_numberOfCellsZ);
}

int CellList::index(const vec3 &position)
{
    int cx = position[0]/m_system->systemSize()[0]*(m_numberOfCellsX-2); // -2 are for the extra ghost cells
    int cy = position[1]/m_system->systemSize()[1]*(m_numberOfCellsY-2);
    int cz = position[2]/m_system->systemSize()[2]*(m_numberOfCellsZ-2);

    return index(cx, cy, cz);
}

void CellList::index3D(const vec3 &position, int &cx, int &cy, int &cz)
{
    cx = position[0]/m_system->systemSize()[0]*(m_numberOfCellsX-2); // -2 are for the extra ghost cells
    cy = position[1]/m_system->systemSize()[1]*(m_numberOfCellsY-2);
    cz = position[2]/m_system->systemSize()[2]*(m_numberOfCellsZ-2);
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

    m_numberOfCellsX = system->systemSize().x() / rCut + 2;
    m_numberOfCellsY = system->systemSize().y() / rCut + 2;
    m_numberOfCellsZ = system->systemSize().z() / rCut + 2;

//    m_numberOfCellsX = system->systemSize().x() / rCut;
//    m_numberOfCellsY = system->systemSize().y() / rCut;
//    m_numberOfCellsZ = system->systemSize().z() / rCut;

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
        Atom *atom = m_system->atoms()[i];
        int cellIndex = index(atom->position);
        m_cells[cellIndex].push_back(atom);
    }

    m_system->resetGhostAtoms();

    for(int cx = 0; cx<=m_numberOfCellsX-2; cx+= m_numberOfCellsX-3) {
        for(int cy = 0; cy<=m_numberOfCellsY-2; cy+= m_numberOfCellsY-3) {
            for(int cz = 0; cz<=m_numberOfCellsZ-2; cz+= m_numberOfCellsZ-3) {
                vector<Atom*> &atoms = m_cells[index(cx,cy,cz)];

                for(int i=0; i<atoms.size(); i++) {
                    Atom *atom = atoms[i];
                    Atom *ghostAtom = m_system->copyAtomToGhostAtom(this, atom);
                    int cellIndex = index(ghostAtom->position);
                    m_cells[cellIndex].push_back(atom);
                }
            }
        }
    }

    CPElapsedTimer::updateCellList().stop();
}
