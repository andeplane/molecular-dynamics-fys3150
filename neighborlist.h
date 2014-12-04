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
    MDDataType_t     m_rShellSquared;
    unsigned int m_numUpdated;
    unsigned int m_totalComparedNeighborPairs;

    unsigned int m_numInitialNeighborPairs;
    unsigned int m_numCurrentNeighborPairs;

    void clear();
public:
    NeighborList();
    void setup(System *system, float rShell);
    void update();
    unsigned int numUpdated() { return m_numUpdated; }
    inline unsigned int *neighborsForAtomWithIndex(int index) { return m_neighbors[index]; }
    CellList &cellList() { return m_cellList; }
    unsigned int totalComparedNeighborPairs() { return m_totalComparedNeighborPairs; }

    MDDataType_t averageNumNeighbors();
    MDDataType_t rShellSquared() { return m_rShellSquared; }

    float currentNeighborPairRatio() { return (float)m_numCurrentNeighborPairs / (float)m_numInitialNeighborPairs; }
    unsigned int &numCurrentNeighborPairs() { return m_numCurrentNeighborPairs; }
};

#endif // NEIGHBORLIST_H
