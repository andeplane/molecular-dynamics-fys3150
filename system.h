#pragma once
#include "vector"
#include "math/vec3.h"
#include "celllist.h"

class Potential; class Integrator;
using std::vector;

class System
{
private:
    vec3 m_systemSize;
    Potential *m_potential;
    Integrator *m_integrator;
    CellList     m_cellList;
    float m_currentTime;
    int m_steps;
    bool m_initialized;
    bool m_shouldSample;
    float m_rCut;
    float m_rShell;
    void validate();
public:
    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, float latticeConstant, float temperature);
    void applyPeriodicBoundaryConditions();
    void removeMomentum();
    void calculateForces();
    void step(float dt);
    unsigned int numberOfAtoms;

    // Setters and getters
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
    void setShouldSample(bool shouldSample);
    void createGhostAtoms();
    void printStatus();
    void setCutoffRadius(float rCut) { m_rCut = rCut; }
    float cutoffRadius() { return m_rCut; }
    CellList &cellList() { return m_cellList; }
};
