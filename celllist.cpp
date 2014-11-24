#include "celllist.h"
#include "system.h"
#include <iostream>
#include "cpelapsedtimer.h"
#include <cassert>
#include <algorithm>
#include <map>
#include <utility>
#include <functional>

using namespace std;

unsigned int CellList::indexf(float x, float y, float z)
{
    unsigned int cx = x*m_oneOverLengthX;
    unsigned int cy = y*m_oneOverLengthY;
    unsigned int cz = z*m_oneOverLengthZ;

    return index(cx, cy, cz);
}

vector<Cell> &CellList::cells()
{
    return m_cells;
}

void CellList::addAtom(vec3 position, vec3 velocity, float mass, unsigned int atomIndex)
{
    unsigned int cellIndex = indexf(position[0], position[1], position[2]);
    assert(cellIndex >= 0 && cellIndex < m_cells.size() && "Trying to add an atom in an out of bounds cell");
    Cell &cell = m_cells[cellIndex];
    cell.addAtom(position, velocity, mass, atomIndex);
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

void CellList::setup(System *system)
{
    m_system = system;
    m_rCut = system->cutoffRadius();

    m_numberOfCellsX = system->systemSize().x() / float(CELLSIZE);
    m_numberOfCellsY = system->systemSize().y() / float(CELLSIZE);
    m_numberOfCellsZ = system->systemSize().z() / float(CELLSIZE);

    m_oneOverLengthX = 1.0 / (system->systemSize().x()/m_numberOfCellsX);
    m_oneOverLengthY = 1.0 / (system->systemSize().y()/m_numberOfCellsY);
    m_oneOverLengthZ = 1.0 / (system->systemSize().z()/m_numberOfCellsZ);

    m_cells.resize(m_numberOfCellsX*m_numberOfCellsY*m_numberOfCellsZ);
    buildNeighborIndicesList();
}


#define make_tuple_with_pair(x, y, z) make_pair(x, make_pair(y,z))
    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

void CellList::buildNeighborIndicesList()
{
    return;
    map<pair<short, pair<short, short> >,bool> cellPairMap;

    short numberOfCellsXHalf = (m_numberOfCellsX-1)/2;
    short numberOfCellsYHalf = (m_numberOfCellsY-1)/2;
    short numberOfCellsZHalf = (m_numberOfCellsZ-1)/2;
    for(short dx = -numberOfCellsXHalf; dx <= numberOfCellsXHalf; dx++) {
        for(short dy = -numberOfCellsYHalf; dy <= numberOfCellsYHalf; dy++) {
            for(short dz = -numberOfCellsZHalf; dz <= numberOfCellsZHalf; dz++) {
                if(cellPairMap.count(make_tuple_with_pair(-dx, -dy, -dz))) continue; // We already have this pair

                // Find the shortest possible distance between two atoms in each cell, if it is larger than rCut, then skip this cell pair
                float dxx = short(dx - sgn(dx))*lengthX();
                float dyy = short(dy - sgn(dy))*lengthY();
                float dzz = short(dz - sgn(dz))*lengthZ();
                float dr2 = dxx*dxx + dyy*dyy + dzz*dzz;
                if(dr2 < m_rCut*m_rCut) {
                    cellPairMap[make_tuple_with_pair(dx,dy,dz)] = true;
                    m_cellNeighborIndicesList.push_back({dx, dy, dz});
                }
            }
        }
    }

    cellPairMap.clear();
}

void CellList::update() {
    CPElapsedTimer::updateCellList().start();
    for(unsigned int cellIndex=0; cellIndex<m_cells.size(); cellIndex++) {
        Cell &oldCell = m_cells[cellIndex];
        for(int atomIndex=0; atomIndex<oldCell.numberOfAtoms; atomIndex++) {
            unsigned int newCellIndex = indexf(oldCell.x[atomIndex], oldCell.y[atomIndex], oldCell.z[atomIndex]);
            if(cellIndex != newCellIndex) {
                // Our atom has moved into another cell
                Cell &newCell = m_cells[newCellIndex];
                // cout << "Atom " << oldCell.index[atomIndex] << " at [" << oldCell.x[atomIndex] << ", " << oldCell.y[atomIndex] << ", " << oldCell.z[atomIndex] << "] moved from cell " << cellIndex << " to " << newCellIndex << " which has " << newCell.numberOfAtoms << " atoms." << endl;
                newCell.addAtomFromCell(oldCell, atomIndex);
                oldCell.removeAtom(atomIndex);
                atomIndex--; // Don't worry if atomIndex is 0, it will wrap around and become UNSIGNED_INT_MAX and then back to zero
            }
        }
    }
    CPElapsedTimer::updateCellList().stop();
}

Cell::Cell() :
    numberOfAtoms(0)
{
    memset(index,0,MAXNUMATOMSPERCELL*sizeof(int));
}

void Cell::addAtomFromCell(Cell &cell, unsigned int i)
{
    addAtom(cell.x[i], cell.y[i], cell.z[i], cell.vx[i], cell.vy[i], cell.vz[i], cell.mass[i], cell.index[i]);
}

void Cell::addAtom(vec3 position, vec3 velocity, float mass, unsigned int atomIndex)
{
    addAtom(position[0], position[1], position[2], velocity[0], velocity[1], velocity[2], mass, atomIndex);
}

void Cell::addAtom(float x, float y, float z, float vx, float vy, float vz, float mass, unsigned int atomIndex)
{
    assert(numberOfAtoms<MAXNUMATOMSPERCELL && "Too many atoms per cell. Increase the MAXNUMATOMSPERCELL size");
    this->x[numberOfAtoms] = x;
    this->y[numberOfAtoms] = y;
    this->z[numberOfAtoms] = z;

    this->vx[numberOfAtoms] = vx;
    this->vy[numberOfAtoms] = vy;
    this->vz[numberOfAtoms] = vz;

    this->mass[numberOfAtoms] = mass;
    this->index[numberOfAtoms] = atomIndex;
    numberOfAtoms++;
}

void Cell::removeAtom(unsigned int i)
{
    if(numberOfAtoms>1) {
        // Switch the atom with the last one and pop the last one
        unsigned int lastAtomIndex = numberOfAtoms-1;
        x[i] = x[lastAtomIndex];
        y[i] = y[lastAtomIndex];
        z[i] = z[lastAtomIndex];
        vx[i] = vx[lastAtomIndex];
        vy[i] = vy[lastAtomIndex];
        vz[i] = vz[lastAtomIndex];
        mass[i] = mass[lastAtomIndex];
        index[i] = index[lastAtomIndex];
    }

    numberOfAtoms--;
}

void CellList::forEachAtom(std::function<void (Cell &cell, unsigned int index)> action)
{
    for(Cell &cell : m_cells) {
        for(unsigned int i=0; i<cell.numberOfAtoms; i++) {
            action(cell, i);
        }
    }
}
