#include "integrators/velocityverlet.h"
#include "system.h"
#include "cpelapsedtimer.h"
#include <iostream>
using namespace std;

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true) // This will set the variable m_firstStep to false when the object is created
{

}

void VelocityVerlet::halfKick(System *system, const float dt)
{
    CPElapsedTimer::halfKick().start();

    system->cellList().forEachAtom([&](Cell &cell, unsigned int i) {
        cell.vx[i] += cell.fx[i]*0.5*dt/cell.mass[i];
        cell.vy[i] += cell.fy[i]*0.5*dt/cell.mass[i];
        cell.vz[i] += cell.fz[i]*0.5*dt/cell.mass[i];
    });
    CPElapsedTimer::halfKick().stop();
}

void VelocityVerlet::move(System *system, const float dt)
{
    CPElapsedTimer::move().start();

    system->cellList().forEachAtom([&](Cell &cell, unsigned int i) {
        cell.x[i] += cell.vx[i]*dt;
        cell.y[i] += cell.vy[i]*dt;
        cell.z[i] += cell.vz[i]*dt;
    });

    CPElapsedTimer::move().stop();
}

void VelocityVerlet::integrate(System *system, float dt)
{
    if(m_firstStep) {
        system->calculateForces();
        m_firstStep = false;
    }
    halfKick(system, dt);
    move(system, dt);
    system->applyPeriodicBoundaryConditions();
    system->calculateForces();
    halfKick(system, dt);
}
