#include "io.h"
#include "system.h"
#include "atom.h"
#include "unitconverter.h"
#include "cstdlib"
#include "cpelapsedtimer.h"
#include <iomanip>
using std::endl; using std::cout;

IO::IO()
{

}

IO::~IO() {
    close();
}

void IO::close() {
    if(m_movieFile.is_open()) m_movieFile.close();
    if(m_statisticsFile.is_open()) m_statisticsFile.close();
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
void IO::saveState(System *system)
{
//    CPElapsedTimer::disk().start();
//    file << system->atoms().numberOfAtoms << endl;
//    file << "The is an optional comment line that can be empty." << endl;
//    for(int n=0; n<system->atoms().numberOfAtoms; n++) {
//        // file << "Ar " << UnitConverter::lengthToAngstroms(atom.position.x()) << " " << UnitConverter::lengthToAngstroms(atom.position.y()) << " " << UnitConverter::lengthToAngstroms(atom.position.z()) << endl;
//    }
//    CPElapsedTimer::disk().stop();
}

void IO::writePerformance(unsigned int timestep) {
    if(!m_performanceFile.is_open()) {
        m_performanceFile.open("performance.txt");

        m_performanceFile << "Timestep     Runtime [s]     DeltaT [s]" << endl;
    }

    m_performanceFile.precision(12);
    m_performanceFile << timestep << "     " << CPElapsedTimer::totalTime() << "      " << CPElapsedTimer::ping() << endl;
}

void IO::writeStatistics(float time, float kineticEnergy, float potentialEnergy, float pressure, float temperature) {
    if(!m_statisticsFile.is_open()) {
        m_statisticsFile.open("statistics.txt");
        m_statisticsFile << "Time [ps]              Kinetic energy [eV]    Potential energy [eV]   Total energy [eV]       Pressure [GPa]         Temperature [K]" << endl;
    }
    m_statisticsFile.precision(12);
    m_statisticsFile << std::scientific << UC::timeToSI(time)*1e12 << "     " << UC::energyToEv(kineticEnergy) << "     " << UC::energyToEv(potentialEnergy) << "     " << UC::energyToEv(kineticEnergy+potentialEnergy) << "     " << UC::pressureToSI(pressure)*1e-9 << "     " << UC::temperatureToSI(temperature) << endl;
}
