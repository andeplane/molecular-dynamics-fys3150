#include <mpi.h>

#include <iostream>
#include <cmath>
#include "math/random.h"

#include "potentials/lennardjones.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include "cpelapsedtimer.h"
#include "modifiers/berendsenthermostat.h"

using namespace std;

void printUnits() {
    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;
}

int main(int args, char *argv[])
{
    MPI_Init(NULL, NULL);
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


    int numTimeSteps = 1e6;
    double dt = UnitConverter::timeFromSI(5e-15); // You should try different values for dt as well.
    int numUnitCells = 12;
    float latticeConstant = 5.26;
    bool loadState = false;
    bool thermostatEnabled = false;
    float temperature = 700;
    // vector<double> temperatures = {100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 555, 560, 565, 570, 575, 580, 585, 590, 595, 600, 605, 610, 615, 620, 625, 630, 635, 640, 650, 660, 670, 680, 690, 700, 800, 900, 1000};
    // vector<double> temperatures = {450, 452, 454, 456, 458, 460, 462, 464, 466, 468, 470, 472, 474, 476, 478, 480, 482, 484, 486, 488, 490, 492, 494, 496, 498, 500, 502, 504, 506, 508, 510, 512};
    vector<double> temperatures = {512, 513, 514, 515, 516, 517, 518, 519, 520, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 533, 534, 535, 536, 537, 538, 539, 540, 541, 542, 543, 544, 545, 546, 547, 548, 549, 550, 551};
    if(world_size > 1) {
        if(world_rank >= temperatures.size()) {
            cout << "Error, too many MPI processes." << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        temperature = temperatures[world_rank];
    }

    if(args>6) {
        dt = UnitConverter::timeFromSI(atof(argv[1])*1e-15);
        numTimeSteps = atoi(argv[2]);
        numUnitCells = atoi(argv[3]);
        latticeConstant = atof(argv[4]);
        loadState = atoi(argv[5]);
        thermostatEnabled = atoi(argv[6]);
        temperature = atof(argv[7]);
    }
    cout << "Starting MD simulation" << endl;
    float rCut = UnitConverter::lengthFromAngstroms(2.5*3.405);

    System system;
    StatisticsSampler statisticsSampler;
    char statisticsFilename[10000];
    sprintf(statisticsFilename, "statistics_%d.txt", int(temperature));
    statisticsSampler.setFilename(string(statisticsFilename));
    BerendsenThermostat thermostat(UnitConverter::temperatureFromSI(temperature), 0.01);

    system.createFCCLattice(numUnitCells, UnitConverter::lengthFromAngstroms(latticeConstant), UnitConverter::temperatureFromSI(temperature));
    system.setPotential(new LennardJones(UnitConverter::lengthFromAngstroms(3.405), 1.0, rCut)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.initialize(rCut);
    system.removeMomentum();

//    IO *movie = new IO(); // To write the state to file
//    movie->open("movie.xyz");

    CPElapsedTimer::timeEvolution().start();

    cout << "Will run " << numTimeSteps << " timesteps." << endl;
    for(int timestep=0; timestep<numTimeSteps; timestep++) {
        bool shouldSample = !(timestep % 100) || thermostatEnabled;
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

        if( !(timestep % 100)) {
            double totalTime = CPElapsedTimer::getInstance().totalTime();
            double timePerTimestep = totalTime / (timestep+1);
            double estimatedTimeLeft = round((numTimeSteps - timestep) * timePerTimestep);
            cout << "Step " << timestep << " Epot/n = " << statisticsSampler.potentialEnergy()/system.atoms().size() << "   Ekin/n = " << statisticsSampler.kineticEnergy()/system.atoms().size() << "   Etot/n = " << statisticsSampler.totalEnergy()/system.atoms().size() << " Time left: " << estimatedTimeLeft << " seconds." << endl;
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
    cout << system.atoms().size()*numTimeSteps / (1000*CPElapsedTimer::totalTime()) << "k atom-timesteps / second. " << endl;

    // movie->close();
    MPI_Finalize();

    return 0;
}
