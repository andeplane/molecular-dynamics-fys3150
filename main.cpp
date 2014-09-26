#include <iostream>
#include <math/random.h>

#include <potentials/lennardjones.h>
#include <integrators/velocityverlet.h>
#include <system.h>
#include <statisticssampler.h>

using namespace std;

int main()
{
    double dt = 0.02;

    System *system = new System();
    system->createFCCLattice(5, 5.26);
    system->setPotential(new LennardJones(1.0, 1.0));
    system->setIntegrator(new VelocityVerlet());

    StatisticsSampler *statisticsSampler = new StatisticsSampler();

    for(int timestep=0; timestep<1000; timestep++) {
        system->step(dt);
        statisticsSampler->sample(system);
    }

    return 0;
}

