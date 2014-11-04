#include <integrators/velocityverlet.h>
#include <system.h>
#include <atom.h>
#include <iostream>
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
    float maxAmplitude = 0;
    for(int i=0; i<system->atoms().size(); i++) {
        Atom *atom = system->atoms()[i];
        atom->velocity += atom->force*(0.5*dt/atom->mass());
        maxAmplitude = std::max(maxAmplitude, atom->force.length());
        if(maxAmplitude > 1e10) {
            cout << "Trouble on atom " << i << endl;
            exit(1);
        }

        // atom->velocity.addAndMultiply(atom->force, 0.5*dt/atom->mass());
    }

    // cout << "Max force: " << maxAmplitude << endl;
    if(maxAmplitude > 1e10) {
        exit(1);
    }
}

void VelocityVerlet::move(System *system, float dt)
{
    float maxAmplitude = 0;
    for(int i=0; i<system->atoms().size(); i++) {
        Atom *atom = system->atoms()[i];
        atom->position += atom->velocity*dt;
        maxAmplitude = std::max(maxAmplitude, atom->velocity.length());
        if(maxAmplitude > 1e10) {
            cout << "Trouble on atom " << i << endl;
            exit(1);
        }
        // atom->position.addAndMultiply(atom->velocity, dt);
    }

    // cout << "Max velocity: " << maxAmplitude << endl;
    if(maxAmplitude > 1e10) {
        exit(1);
    }
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
