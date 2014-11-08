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
    vector<vector<int> > m_cells;

public:
    CellList();
    void setup(System *system, float rCut);
    void clear();
    void update();
    int numberOfCellsX() { return m_numberOfCellsX; }
    int numberOfCellsY() { return m_numberOfCellsY; }
    int numberOfCellsZ() { return m_numberOfCellsZ; }
    int index(int cx, int cy, int cz);
    int indexPeriodic(int cx, int cy, int cz);
    int index(const vec3 &position);
    vector<int> &operator[](int index) { return m_cells[index]; }
    vector<vector<int> > &cells();
};

#endif // CELLLIST_H
