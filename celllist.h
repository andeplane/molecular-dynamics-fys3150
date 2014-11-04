#ifndef CELLLIST_H
#define CELLLIST_H
#include <vector>
using std::vector;
#include "math/vec3.h"
#include "system.h"

class Atom; class System;
class CellList
{
private:
    System *m_system;
    int m_numberOfCellsX;
    int m_numberOfCellsY;
    int m_numberOfCellsZ;
    vector<vector<Atom*> > m_cells;
    inline int index(int cx, int cy, int cz) {
        return cx*m_numberOfCellsY*m_numberOfCellsZ + cy*m_numberOfCellsZ + cz;
    }

    inline int index(const vec3 &position) {
        int cx = position[0]/m_system->systemSize()[0];
        int cy = position[1]/m_system->systemSize()[1];
        int cz = position[2]/m_system->systemSize()[2];

        return index(cx, cy, cz);
    }

public:
    CellList();
    void setup(System *system, float rCut);
    void clear();
    void update();
};

#endif // CELLLIST_H
