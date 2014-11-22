#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H
class System;
#include "celllist.h"
#include <vector>
using std::vector;

class NeighborList
{
private:
    vector<vector<Atom*> > m_neighbors;
    CellList  m_cellList;
    System   *m_system;
    float     m_rShellSquared;

    void clear();
public:
    NeighborList();
    void setup(System *system, float rShell);
    void update();
    vector<Atom*> &neighborsForAtomWithIndex(int index) { return m_neighbors[index]; }
};

#endif // NEIGHBORLIST_H
