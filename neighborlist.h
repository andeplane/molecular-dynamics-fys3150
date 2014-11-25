#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H
class System;
#include "celllist.h"
#include "atoms.h"
#include <vector>
using std::vector;

class NeighborList
{
private:
    // vector<vector<unsigned int> > m_neighbors;
    unsigned int **m_neighbors;
    CellList  m_cellList;
    System   *m_system;
    float     m_rShellSquared;
    unsigned int m_numNeighborPairs;

    void clear();
public:
    NeighborList();
    void setup(System *system, float rShell);
    void update();
    inline unsigned int *neighborsForAtomWithIndex(int index) { return m_neighbors[index]; }
    CellList &cellList() { return m_cellList; }
    unsigned int numNeighborPairs() { return m_numNeighborPairs; }
    float averageNumNeighbors();
};

#endif // NEIGHBORLIST_H
