#include "celllist.h"
#include "system.h"
#include <iostream>
#include "cpelapsedtimer.h"
#include <cassert>
#include <algorithm>

using namespace std;

unsigned int CellList::indexf(float x, float y, float z)
{
    unsigned int cx = double(x)*m_oneOverLengthX;
    unsigned int cy = double(y)*m_oneOverLengthY;
    unsigned int cz = double(z)*m_oneOverLengthZ;

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

void CellList::setup(System *system, float rCut)
{
    m_system = system;

    m_numberOfCellsX = system->systemSize().x() / rCut;
    m_numberOfCellsY = system->systemSize().y() / rCut;
    m_numberOfCellsZ = system->systemSize().z() / rCut;

    m_oneOverLengthX = 1.0 / (system->systemSize().x()/m_numberOfCellsX);
    m_oneOverLengthY = 1.0 / (system->systemSize().y()/m_numberOfCellsY);
    m_oneOverLengthZ = 1.0 / (system->systemSize().z()/m_numberOfCellsZ);

    m_cells.resize(m_numberOfCellsX*m_numberOfCellsY*m_numberOfCellsZ);
}

void CellList::update() {
    for(unsigned int cellIndex=0; cellIndex<m_cells.size(); cellIndex++) {
        Cell &oldCell = m_cells[cellIndex];
        for(unsigned int atomIndex=0; atomIndex<oldCell.numberOfAtoms; atomIndex++) {
            unsigned int newCellIndex = indexf(oldCell.x[atomIndex], oldCell.y[atomIndex], oldCell.z[atomIndex]);
            if(cellIndex != newCellIndex) {
                // Our atom has moved into another cell
                Cell &newCell = m_cells[newCellIndex];
                newCell.addAtomFromCell(oldCell, atomIndex);
                oldCell.removeAtom(atomIndex);
                atomIndex--;
            }
        }
    }
}

Cell::Cell() :
    numberOfAtoms(0)
{
    memset(cellIndex,0,MAXNUMATOMS*sizeof(int));
    memset(index,0,MAXNUMATOMS*sizeof(int));
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
        unsigned int lastAtomIndex = numberOfAtoms-1;
        x[i] = x[lastAtomIndex];
        y[i] = y[lastAtomIndex];
        z[i] = z[lastAtomIndex];
        vx[i] = vx[lastAtomIndex];
        vy[i] = vy[lastAtomIndex];
        vz[i] = vz[lastAtomIndex];
        mass[i] = mass[lastAtomIndex];
        index[i] = index[lastAtomIndex];
        numberOfAtoms--;
    }
}

void CellList::forEachAtom(std::function<void (Cell &cell, unsigned int index)> action)
{
    for(Cell &cell : m_cells) {
        for(unsigned int i=0; i<cell.numberOfAtoms; i++) {
            action(cell, i);
        }
    }
}
