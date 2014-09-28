#include <iostream>
#include <math/random.h>

#include <potentials/lennardjones.h>
#include <integrators/eulercromer.h>
#include <system.h>
#include <statisticssampler.h>
#include <atom.h>
#include <io.h>
#include <unitconverter.h>

using namespace std;

int main()
{
    double dt = 100;

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;

    System *system = new System();
    system->createFCCLattice(5, 5.26);
    system->setPotential(new LennardJones(1.0, 1.0));
    system->setIntegrator(new EulerCromer());
    system->removeMomentum();

    for(int n=0; n<100; n++) {
        // Add one example atom. You'll have to create many such atoms in the createFCCLattice function above.
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26)); // Argon mass, see http://en.wikipedia.org/wiki/Argon
        atom->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
        atom->position.randomUniform(UnitConverter::lengthFromAngstroms(-10), UnitConverter::lengthFromAngstroms(10));
        system->atoms().push_back(atom); // Add it to the list of atoms
    }

    StatisticsSampler *statisticsSampler = new StatisticsSampler(); //

    IO *movie = new IO(); // To write the state to file
    movie->open("movie.xyz");

    for(int timestep=0; timestep<1000; timestep++) {
        system->step(dt);
        statisticsSampler->sample(system);

        movie->saveState(system);
    }

    movie->close();

    return 0;
}
