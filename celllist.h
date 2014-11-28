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
    float m_oneOverLengthX;
    float m_oneOverLengthY;
    float m_oneOverLengthZ;
    vector<vector<unsigned int> > m_cells;
    vector<vector<vector<unsigned int>*> > m_neighbors;

public:
    CellList();
    void setup(System *system, float rCut);
    void clear();
    void update();
    float lengthX() { return 1.0/m_oneOverLengthX; }
    float lengthY() { return 1.0/m_oneOverLengthY; }
    float lengthZ() { return 1.0/m_oneOverLengthZ; }

    int numberOfCellsX() { return m_numberOfCellsX; }
    int numberOfCellsY() { return m_numberOfCellsY; }
    int numberOfCellsZ() { return m_numberOfCellsZ; }

    void index3d(unsigned int &index, unsigned int &i, unsigned int &j, unsigned int &k) {
        k = (index / (m_numberOfCellsX*m_numberOfCellsY));
        j = (index / m_numberOfCellsX) % m_numberOfCellsY;
        i = index % m_numberOfCellsX;
    }

    unsigned int indexf(float x, float y, float z);
    unsigned int index(unsigned int cx, unsigned int cy, unsigned int cz) { return cx + cy*m_numberOfCellsX + cz*m_numberOfCellsY*m_numberOfCellsX; }
    unsigned int indexPeriodic(int cx, int cy, int cz) { return ( (cx+m_numberOfCellsX) % m_numberOfCellsX) + ( (cy+m_numberOfCellsY) % m_numberOfCellsY)*m_numberOfCellsX + ( (cz+m_numberOfCellsZ) % m_numberOfCellsZ)*m_numberOfCellsX*m_numberOfCellsY; }
    vector<vector<unsigned int> *> &getNeighbors(unsigned int cx, unsigned int cy, unsigned int cz) { return m_neighbors[index(cx,cy,cz)]; }
    vector<vector<unsigned int> *> &getNeighbors(unsigned int cellIndex) { return m_neighbors[cellIndex]; }
    vector<unsigned int> &operator[](int index) { return m_cells[index]; }
    vector<vector<unsigned int> > &cells();
};

#endif // CELLLIST_H
