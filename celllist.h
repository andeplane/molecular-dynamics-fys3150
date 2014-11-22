#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
using std::vector;
#include "math/vec3.h"

class Atom; class System;
class CellList
{
private:
    System *m_system;
    int m_numberOfCellsX;
    int m_numberOfCellsY;
    int m_numberOfCellsZ;
    vector<vector<Atom*> > m_cells;

public:
    CellList();
    void setup(System *system, float rCut);
    void clear();
    void update();
    int numberOfCellsX() { return m_numberOfCellsX; }
    int numberOfCellsY() { return m_numberOfCellsY; }
    int numberOfCellsZ() { return m_numberOfCellsZ; }
    int index(int cx, int cy, int cz) { return cx*m_numberOfCellsY*m_numberOfCellsZ + cy*m_numberOfCellsZ + cz; }
    int indexPeriodic(int cx, int cy, int cz) { return ( (cx+m_numberOfCellsX) % m_numberOfCellsX)*m_numberOfCellsY*m_numberOfCellsZ + ( (cy+m_numberOfCellsY) % m_numberOfCellsY)*m_numberOfCellsZ + ( (cz+m_numberOfCellsZ) % m_numberOfCellsZ); }
    int index(const vec3 &position);
    void sort();
    vector<Atom*> &operator[](int index) { return m_cells[index]; }
    vector<vector<Atom*> > &cells();
};

#endif // CELLLIST_H
