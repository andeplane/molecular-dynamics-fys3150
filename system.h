#pragma once
#include "vector"
#include "math/vec3.h"
#include "celllist.h"
#include "neighborlist.h"
#include "atoms.h"

class Potential; class Integrator;
using std::vector;

class System
{
private:
    vec3 m_systemSize;
    Atoms *m_atoms;
    MiniAtoms *m_miniAtoms;
    Potential *m_potential;
    Integrator *m_integrator;

    NeighborList m_neighborList;
    double m_currentTime;
    int m_steps;
    bool m_initialized;
    bool m_shouldSample;
    MDDataType_t m_rCut;
    MDDataType_t m_rShell;
    void validate();
public:
    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, float latticeConstant, float temperature);
    void applyPeriodicBoundaryConditions();
    void removeMomentum();
    void calculateForces();
    void step(double dt);

    // Setters and getters
    Atoms &atoms() { return *m_atoms; }
    MiniAtoms &miniAtoms() { return *m_miniAtoms; }
    vec3 &systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double currentTime() { return m_currentTime; }
    void setCurrentTime(float currentTime) { m_currentTime = currentTime; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
    double volume() { return m_systemSize[0]*m_systemSize[1]*m_systemSize[2]; }
    NeighborList &neighborList() { return m_neighborList; }
    void initialize(MDDataType_t cutoffRadius);
    void setShouldSample(bool shouldSample);
    void createGhostAtoms();
    void printStatus();
    MDDataType_t rCut() { return m_rCut; }
    bool initialized() const;
};
