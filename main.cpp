#include <iostream>
#include <math/random.h>

#include <potentials/lennardjones.h>
#include <integrators/velocityverlet.h>
#include <system.h>
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include "cpelapsedtimer.h"

using namespace std;

int main()
{
    int numTimeSteps = 20;
    double dt = UnitConverter::timeFromSI(1e-14); // You should try different values for dt as well.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;

    float rCut = UnitConverter::lengthFromAngstroms(2.5*3.405);
    System system;
    system.createFCCLattice(10, UnitConverter::lengthFromAngstroms(5.26*1.0875335), UnitConverter::temperatureFromSI(300));
    system.setPotential(new LennardJones(UnitConverter::lengthFromAngstroms(3.405), 1.0, rCut)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.initialize(rCut);
    system.removeMomentum();

    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //

    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    for(int timestep=0; timestep<numTimeSteps; timestep++) {
        system.step(dt);
        if( !(timestep % 2)) {
            // statisticsSampler->sample(&system);
            // cout << "Step " << timestep << " Epot = " << statisticsSampler->samplePotentialEnergy(&system) << "   Ekin = " << statisticsSampler->sampleKineticEnergy(&system) << "   Etot = " << statisticsSampler->totalEnergy() <<  endl;
            cout << "Step " << timestep << endl;
        }
        // movie->saveState(&system);
    }

    float calculateForcesFraction = CPElapsedTimer::calculateForces().elapsedTime() / CPElapsedTimer::totalTime();
    float halfKickFraction = CPElapsedTimer::halfKick().elapsedTime() / CPElapsedTimer::totalTime();
    float moveFraction = CPElapsedTimer::move().elapsedTime() / CPElapsedTimer::totalTime();
    float updateNeighborListFraction = CPElapsedTimer::updateNeighborList().elapsedTime() / CPElapsedTimer::totalTime();
    float updateCellListFraction = CPElapsedTimer::updateCellList().elapsedTime() / CPElapsedTimer::totalTime();
    float periodicBoundaryConditionsFraction = CPElapsedTimer::periodicBoundaryConditions().elapsedTime() / CPElapsedTimer::totalTime();
    float samplingFraction = CPElapsedTimer::sampling().elapsedTime() / CPElapsedTimer::totalTime();

    cout << endl << "Program finished after " << CPElapsedTimer::totalTime() << " seconds. Time analysis:" << endl;
    cout << fixed
         << "      Force calculation : " << CPElapsedTimer::calculateForces().elapsedTime() << " s ( " << 100*calculateForcesFraction << "%)" <<  endl
         << "      Moving            : " << CPElapsedTimer::move().elapsedTime() << " s ( " << 100*moveFraction << "%)" <<  endl
         << "      Half kick         : " << CPElapsedTimer::halfKick().elapsedTime() << " s ( " << 100*halfKickFraction << "%)" <<  endl
         << "      Update neighbors  : " << CPElapsedTimer::updateNeighborList().elapsedTime() << " s ( " << 100*updateNeighborListFraction << "%)" <<  endl
         << "      Update cells      : " << CPElapsedTimer::updateCellList().elapsedTime() << " s ( " << 100*updateCellListFraction << "%)" <<  endl
         << "      Periodic boundary : " << CPElapsedTimer::periodicBoundaryConditions().elapsedTime() << " s ( " << 100*periodicBoundaryConditionsFraction << "%)" <<  endl
         << "      Sampling          : " << CPElapsedTimer::sampling().elapsedTime() << " s ( " << 100*samplingFraction << "%)" <<  endl;
    cout << endl << numTimeSteps / CPElapsedTimer::totalTime() << " timesteps / second. " << endl;
    cout << system.atoms().size()*numTimeSteps / (1000*CPElapsedTimer::totalTime()) << "k atom-timesteps / second. " << endl;

    movie->close();

    return 0;
}
