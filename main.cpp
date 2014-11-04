#include <iostream>
#include <math/random.h>

#include <potentials/lennardjones.h>
#include <integrators/velocityverlet.h>
#include <system.h>
#include <statisticssampler.h>
#include <atom.h>
#include <io.h>
#include <unitconverter.h>

using namespace std;

void testVec3() {
    vec3 test(3,2,1);
    cout << "test = " << test << endl;
    test += 2;
    cout << "+= 2: " << test << endl;

    test -= 2;
    cout << "-= 2: " << test << endl;

    test *= 2;
    cout << "*= 2: " << test << endl;

    test /= 2;
    cout << "/= 2: " << test << endl;

    test += vec3(2,5,10);
    cout << "+= vec3(2,5,10); " << test << endl;

    test -= vec3(2,5,10);
    cout << "-= vec3(2,5,10); " << test << endl;

    test *= vec3(2,5,10);
    cout << "*= vec3(2,5,10); " << test << endl;

    test /= vec3(2,5,10);
    cout << "/= vec3(2,5,10); " << test << endl;
    vec3 test2 = -test;
    cout << "-test; " << test2 << endl;

    exit(1);
}

int main()
{

    double dt = UnitConverter::timeFromSI(1e-14); // You should try different values for dt as well.

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
    cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;

    float rCut = UnitConverter::lengthFromAngstroms(2.5*3.405);
    System system;
    system.createFCCLattice(10, UnitConverter::lengthFromAngstroms(5.26), UnitConverter::temperatureFromSI(300));
    system.setPotential(new LennardJones(UnitConverter::lengthFromAngstroms(3.405), 1.0, rCut)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.initialize(rCut);
    system.removeMomentum();

    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //

    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    for(int timestep=0; timestep<1000; timestep++) {
        system.step(dt);
        if( !(timestep % 100)) {
            // statisticsSampler->sample(&system);

            cout << "Step " << timestep << " Epot = " << statisticsSampler->samplePotentialEnergy(&system) << "   Ekin = " << statisticsSampler->sampleKineticEnergy(&system) << "   Etot = " << statisticsSampler->totalEnergy() <<  endl;
        }
        // movie->saveState(&system);
    }

    movie->close();

    return 0;
}
