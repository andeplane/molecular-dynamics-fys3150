#include <integrators/velocityverlet.h>
#include <system.h>
#include <iostream>
#include "cpelapsedtimer.h"
using namespace std;

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true) // This will set the variable m_firstStep to false when the object is created
{

}

void VelocityVerlet::halfKick(System *system, const float dt)
{
    CPElapsedTimer::halfKick().start();
    Atoms &atoms = system->atoms();
    const float dtHalf = dt*0.5;
#ifdef MD_SIMD
#pragma simd
#endif
    for(int i=0; i<system->atoms().numberOfAtoms; i++) {
        atoms.vx[i] += atoms.fx[i]*dtHalf*atoms.inverseMass[i];
        atoms.vy[i] += atoms.fy[i]*dtHalf*atoms.inverseMass[i];
        atoms.vz[i] += atoms.fz[i]*dtHalf*atoms.inverseMass[i];
    }
    CPElapsedTimer::halfKick().stop();
}

void VelocityVerlet::move(System *system, const float dt)
{
    CPElapsedTimer::move().start();
    Atoms &atoms = system->atoms();
#ifdef MD_SIMD
#pragma simd
#endif
    for(int i=0; i<system->atoms().numberOfAtoms; i++) {
        atoms.x[i] += atoms.vx[i]*dt;
        atoms.y[i] += atoms.vy[i]*dt;
        atoms.z[i] += atoms.vz[i]*dt;
    }
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
