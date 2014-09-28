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
    double dt = 0.02;

    System *system = new System();
    system->createFCCLattice(5, 5.26);
    system->setPotential(new LennardJones(1.0, 1.0));
    system->setIntegrator(new EulerCromer());
    system->removeMomentum();

    // Add one example atom. You'll have to create many such atoms in the createFCCLattice function above.
    Atom *atom = new Atom(39.948); // Argon mass in atomic units, see http://en.wikipedia.org/wiki/Argon
    atom->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(300));
    system->atoms().push_back(atom); // Add it to the list of atoms

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
