#include <iostream>
#include <cstdlib>
#define ALIGNMALLOC 64
#include "math/random.h"

#include "potentials/lennardjones.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "io.h"
#include "unitconverter.h"
#include "cpelapsedtimer.h"
#include "modifiers/berendsenthermostat.h"

using namespace std;

unsigned long calculateFlops(System *system, unsigned int numTimesteps) {
    unsigned long flops = system->atoms().numberOfComputedForces*45 + system->neighborList().totalComparedNeighborPairs()*22; // Forces and neighbor list builds
    flops += system->atoms().numberOfAtoms*numTimesteps*(9 + 6 + 9); // Velocity verlet
    flops += system->atoms().numberOfAtoms*numTimesteps*(18); // Periodic boundary conditions
    return flops;
}

System mdSystem;
IO fileHandler; // To write the state to file
StatisticsSampler statisticsSampler(&fileHandler);
double dt = UnitConverter::timeFromSI(5e-15);
unsigned int timestep = 0;
unsigned int measureEvery = 1;
extern "C" {

double *x() {
    return mdSystem.atoms().x;
}

double *y() {
    return mdSystem.atoms().y;
}

double *z() {
    return mdSystem.atoms().z;
}

int numberOfAtoms() {
    return mdSystem.atoms().numberOfAtoms;
}

double *systemSize() {
    return &mdSystem.systemSize()[0];
}

void setMeasureEvery(unsigned int val) {
    measureEvery = val;
}

void initialize(int numUnitCells = 8, float temperature = 150) {
    float latticeConstant = 5.26;
    float rCut = UnitConverter::lengthFromAngstroms(2.5*3.405);

    mdSystem.createFCCLattice(numUnitCells, UnitConverter::lengthFromAngstroms(latticeConstant), UnitConverter::temperatureFromSI(temperature));
    mdSystem.setPotential(new LennardJones(UnitConverter::lengthFromAngstroms(3.405), 1.0, rCut)); // You must insert correct parameters here
    mdSystem.setIntegrator(new VelocityVerlet());
    mdSystem.initialize(rCut);
    mdSystem.removeMomentum();
}

void step(int timesteps = 1) {
    if(!mdSystem.initialized()) {
        initialize();
    }

    for(unsigned int i=0; i<timesteps; i++) {
        bool shouldSample = !(timestep % measureEvery);
        mdSystem.setShouldSample(shouldSample);
        mdSystem.step(dt);

        if(shouldSample) {
            statisticsSampler.sample(&mdSystem);
        }

        if( shouldSample) {
            cout << "Step " << timestep << " t= " << UnitConverter::timeToSI(mdSystem.currentTime())*1e12 << " ps   Epot/n = " << statisticsSampler.potentialEnergy()/mdSystem.atoms().numberOfAtoms << "   Ekin/n = " << statisticsSampler.kineticEnergy()/mdSystem.atoms().numberOfAtoms << "   Etot/n = " << statisticsSampler.totalEnergy()/mdSystem.atoms().numberOfAtoms <<  endl;
        }

        timestep++;
    }
}
}

int main(int args, char *argv[])
{

    return 0;
}
