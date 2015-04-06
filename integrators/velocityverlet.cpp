#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"
#include <iostream>
#include "../cpelapsedtimer.h"
using namespace std;

VelocityVerlet::VelocityVerlet() :
    m_firstStep(true) // This will set the variable m_firstStep to false when the object is created
{

}

VelocityVerlet::~VelocityVerlet()
{

}

void VelocityVerlet::halfKick(System *system, float dt)
{
    CPElapsedTimer::halfKick().start();
    #pragma ivdep
    for(int i=0; i<system->atoms().size(); i++) {
        Atom &atom = system->atoms()[i];
        atom.velocity.addAndMultiply(atom.force, 0.5*dt/atom.mass());
    }
    CPElapsedTimer::halfKick().stop();
}

void VelocityVerlet::move(System *system, float dt)
{
    CPElapsedTimer::move().start();
    #pragma ivdep
    for(int i=0; i<system->atoms().size(); i++) {
        Atom &atom = system->atoms()[i];
        atom.position.addAndMultiply(atom.velocity, dt);
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
