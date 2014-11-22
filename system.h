#pragma once
#include "vector"
#include "atom.h"
#include "math/vec3.h"
#include "celllist.h"
#include "neighborlist.h"

class Potential; class Integrator;
using std::vector;

class System
{
private:
    vec3 m_systemSize;
    vector<Atom> m_atoms;
    Potential *m_potential;
    Integrator *m_integrator;

    CellList m_cellList;
    NeighborList m_neighborList;
    float m_currentTime;
    int m_steps;
    bool m_initialized;
public:
    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, float latticeConstant, float temperature);
    void applyPeriodicBoundaryConditions();
    void removeMomentum();
    void calculateForces();
    void step(float dt);

    // Setters and getters
    vector<Atom> &atoms() { return m_atoms; }
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    float currentTime() { return m_currentTime; }
    void setCurrentTime(float currentTime) { m_currentTime = currentTime; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    float volume() { return m_systemSize[0]*m_systemSize[1]*m_systemSize[2]; }
    CellList &cellList() { return m_cellList; }
    NeighborList &neighborList() { return m_neighborList; }
    void initialize(float cutoffRadius);
    void setShouldSample(bool shouldSample);
};
