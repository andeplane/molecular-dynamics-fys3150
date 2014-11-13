#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H
class System; class Atom;
#include "celllist.h"
#include <vector>
using std::vector;

class NeighborList
{
private:
    vector<vector<float*> > m_neighbors;
    CellList  m_cellList;
    System   *m_system;
    float     m_rShellSquared;
    vec3      m_deltaRVector;

    void clear();
public:
    NeighborList();
    void setup(System *system, float rShell);
    void update();
    vector<float*> &neighborsForAtomWithIndex(int index) { return m_neighbors[index]; }
};

#endif // NEIGHBORLIST_H
