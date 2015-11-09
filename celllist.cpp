#include "celllist.h"
#include "system.h"
#include "atom.h"
#include <iostream>
using namespace std;

CellList::CellList()
{

}

void CellList::build(System *system)
{
    cout << "Cell list building is disabled because System::applyPeriodicBoundaryConditions is not implemented yet" << endl;
    return;
    int nx = system->systemSize().x() / m_cutoffDistance;
    int ny = system->systemSize().y() / m_cutoffDistance;
    int nz = system->systemSize().z() / m_cutoffDistance;
    m_cutoffDistance = system->systemSize().x() / nx;
    m_cutoffDistance = std::max(m_cutoffDistance, system->systemSize().y() / ny);
    m_cutoffDistance = std::max(m_cutoffDistance, system->systemSize().z() / nz);
    // Now with the new (correct) cutoff distance, calculate the number of cells.
    nx = system->systemSize().x() / m_cutoffDistance;
    ny = system->systemSize().y() / m_cutoffDistance;
    nz = system->systemSize().z() / m_cutoffDistance;

    // Resize so our cell vectors have the correct size
    m_cells.resize(nx);
    for(int i=0; i<nx; i++) {
        m_cells[i].resize(ny);
        for(int j=0; j<ny; j++) {
            m_cells[i][j].resize(nz);
            for(int k=0; k<nz; k++) {
                m_cells[i][j][k].clear();
            }
        }
    }

    // Compute the cell dimensions
    double cellSizeX = system->systemSize().x()/nx;
    double cellSizeY = system->systemSize().y()/ny;
    double cellSizeZ = system->systemSize().z()/nz;
    for(Atom *atom : system->atoms()) {
        int i = atom->position.x()/cellSizeX;
        int j = atom->position.y()/cellSizeY;
        int k = atom->position.z()/cellSizeZ;
        m_cells.at(i).at(j).at(k).push_back(atom);
    }
}

