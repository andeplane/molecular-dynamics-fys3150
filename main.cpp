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
    unsigned long flops = system->atoms().numberOfComputedForces*45 + system->neighborList().numNeighborPairs()*22; // Forces and neighbor list builds
    flops += system->atoms().numberOfAtoms*numTimesteps*(9 + 6 + 9); // Velocity verlet
    flops += system->atoms().numberOfAtoms*numTimesteps*(18); // Periodic boundary conditions
    return flops;
}

int main(int args, char *argv[])
{
    unsigned int numTimeSteps = 1e8;
    double dt = UnitConverter::timeFromSI(1e-14); // You should try different values for dt as well.
    int numUnitCells = 6;
    float latticeConstant = 5.26;
    // float latticeConstant = 5.885;
    bool loadState = false;
    bool thermostatEnabled = false;
    float temperature = 1000;
    if(args>1) {
        dt = UnitConverter::timeFromSI(atof(argv[1])*1e-15);
        numTimeSteps = atoi(argv[2]);
        numUnitCells = atoi(argv[3]);
        latticeConstant = atof(argv[4]);
        loadState = atoi(argv[5]);
        thermostatEnabled = atoi(argv[6]);
        temperature = atof(argv[7]);
    }

    float rCut = UnitConverter::lengthFromAngstroms(2.5*3.405);

    System system;
    IO fileHandler; // To write the state to file
    StatisticsSampler statisticsSampler(&fileHandler);
    BerendsenThermostat thermostat(UnitConverter::temperatureFromSI(temperature), 0.01);

    system.createFCCLattice(numUnitCells, UnitConverter::lengthFromAngstroms(latticeConstant), UnitConverter::temperatureFromSI(temperature));
    system.setPotential(new LennardJones(UnitConverter::lengthFromAngstroms(3.405), 1.0, rCut)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.initialize(rCut);
    system.removeMomentum();

    int measureEvery = 1000;

    CPElapsedTimer::timeEvolution().start();
    cout << "Will run " << numTimeSteps << " timesteps." << endl;
    for(int timestep=0; timestep<numTimeSteps; timestep++) {
        bool shouldSample = !(timestep % measureEvery) || thermostatEnabled;
        system.setShouldSample(shouldSample);
        system.step(dt);

        if(shouldSample) {
            CPElapsedTimer::sampling().start();
            statisticsSampler.sample(&system);
            CPElapsedTimer::sampling().stop();
        }

        if(thermostatEnabled) {
            CPElapsedTimer::thermostat().start();
            thermostat.apply(&system, &statisticsSampler);
            CPElapsedTimer::thermostat().stop();
        }

        if( !(timestep % measureEvery)) {
            cout << "Step " << timestep << " t= " << UnitConverter::timeToSI(system.currentTime())*1e12 << " ps   Epot/n = " << statisticsSampler.potentialEnergy()/system.atoms().numberOfAtoms << "   Ekin/n = " << statisticsSampler.kineticEnergy()/system.atoms().numberOfAtoms << "   Etot/n = " << statisticsSampler.totalEnergy()/system.atoms().numberOfAtoms <<  endl;
            fileHandler.writePerformance(timestep);
        }
        // movie->saveState(&system);
    }
    CPElapsedTimer::timeEvolution().stop();


    float calculateForcesFraction = CPElapsedTimer::calculateForces().elapsedTime() / CPElapsedTimer::totalTime();
    float halfKickFraction = CPElapsedTimer::halfKick().elapsedTime() / CPElapsedTimer::totalTime();
    float moveFraction = CPElapsedTimer::move().elapsedTime() / CPElapsedTimer::totalTime();
    float updateNeighborListFraction = CPElapsedTimer::updateNeighborList().elapsedTime() / CPElapsedTimer::totalTime();
    float updateCellListFraction = CPElapsedTimer::updateCellList().elapsedTime() / CPElapsedTimer::totalTime();
    float periodicBoundaryConditionsFraction = CPElapsedTimer::periodicBoundaryConditions().elapsedTime() / CPElapsedTimer::totalTime();
    float samplingFraction = CPElapsedTimer::sampling().elapsedTime() / CPElapsedTimer::totalTime();
    float timeEvolutionFraction = CPElapsedTimer::timeEvolution().elapsedTime() / CPElapsedTimer::totalTime();
    float thermostatFraction = CPElapsedTimer::thermostat().elapsedTime() / CPElapsedTimer::totalTime();

    cout << endl << "Program finished after " << CPElapsedTimer::totalTime() << " seconds. Time analysis:" << endl;
    cout << fixed
         << "      Time evolution    : " << CPElapsedTimer::timeEvolution().elapsedTime() << " s ( " << 100*timeEvolutionFraction << "%)" <<  endl
         << "      Force calculation : " << CPElapsedTimer::calculateForces().elapsedTime() << " s ( " << 100*calculateForcesFraction << "%)" <<  endl
         << "      Thermostat        : " << CPElapsedTimer::thermostat().elapsedTime() << " s ( " << 100*thermostatFraction << "%)" <<  endl
         << "      Moving            : " << CPElapsedTimer::move().elapsedTime() << " s ( " << 100*moveFraction << "%)" <<  endl
         << "      Half kick         : " << CPElapsedTimer::halfKick().elapsedTime() << " s ( " << 100*halfKickFraction << "%)" <<  endl
         << "      Update neighbors  : " << CPElapsedTimer::updateNeighborList().elapsedTime() << " s ( " << 100*updateNeighborListFraction << "%)" <<  endl
         << "      Update cells      : " << CPElapsedTimer::updateCellList().elapsedTime() << " s ( " << 100*updateCellListFraction << "%)" <<  endl
         << "      Periodic boundary : " << CPElapsedTimer::periodicBoundaryConditions().elapsedTime() << " s ( " << 100*periodicBoundaryConditionsFraction << "%)" <<  endl
         << "      Sampling          : " << CPElapsedTimer::sampling().elapsedTime() << " s ( " << 100*samplingFraction << "%)" <<  endl;
    cout << endl << numTimeSteps / CPElapsedTimer::totalTime() << " timesteps / second. " << endl;
    cout << system.atoms().numberOfAtoms*numTimeSteps / (1000*CPElapsedTimer::totalTime()) << "k atom-timesteps / second. " << endl;
    cout << "Average number of neighbors per atom: " << system.neighborList().averageNumNeighbors() << endl;
    float totalTimePerDay = dt*numTimeSteps/CPElapsedTimer::totalTime() * 86400;
    float nanoSecondsPerDay = UnitConverter::timeToSI(totalTimePerDay)*1e9;
    cout << "Estimated " << nanoSecondsPerDay << " ns simulated time per day" << endl;

    // Performance numbers
    int pairsPerSecond = system.atoms().numberOfComputedForces / CPElapsedTimer::totalTime();
    int neighborPairsPerSecond = system.neighborList().numNeighborPairs() / CPElapsedTimer::totalTime();
    unsigned long flops = calculateFlops(&system, numTimeSteps);
    float flopsPerSecond = flops / CPElapsedTimer::totalTime();
    unsigned long bytesPerSecond = (pairsPerSecond*12 + neighborPairsPerSecond*6)*sizeof(float);
    cout << pairsPerSecond/1e6 << " mega pairs computed per second (" << system.atoms().numberOfComputedForces/1e6 << " mega pairs total)" << endl;
    cout << neighborPairsPerSecond/1e6 << " mega neighbor pairs computed per second (" << neighborPairsPerSecond/1e6 << " mega pairs total)" << endl;
    cout << "Memory read speed: " << bytesPerSecond/1e9 << " gigabytes / sec." << endl;
    cout << "Flops: " << flopsPerSecond/1e9 << " Gflops / sec." << endl;

    return 0;
}
