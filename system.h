#pragma once
#include <vector>
#include <atom.h>
#include <math/vec3.h>

class Potential; class Integrator;
using std::vector;
using CompPhys::vec3;

class System
{
private:
    vec3 m_systemSize;
    vector<Atom*> m_atoms;
    Potential *m_potential;
    Integrator *m_integrator;
    double m_currentTime;
    int m_steps;

public:
    System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCells, double latticeConstant);
    void applyPeriodicBoundaryConditions();
    void step(double dt);

    // Setters and getters
    vector<Atom *> &atoms() { return m_atoms; }
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 &systemSize) { m_systemSize = systemSize; }
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double currentTime() { return m_currentTime; }
    void setCurrentTime(double currentTime) { m_currentTime = currentTime; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
};
