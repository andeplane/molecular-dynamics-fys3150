#include "io.h"
#include "system.h"
#include "unitconverter.h"
#include "cstdlib"
#include "cpelapsedtimer.h"
#include <iomanip>
#include <fstream>
using std::endl; using std::cout; using std::ofstream;

IO::IO()
{
    m_stateFileIndex = 0;
    filename = new char[1000];
}

IO::~IO() {
    close();
}

void IO::close() {
    if(m_movieFile.is_open()) m_movieFile.close();
    if(m_statisticsFile.is_open()) m_statisticsFile.close();
}

void IO::savePositionsBinary(System *system) {
    sprintf(filename, "states/%04d.bin", m_stateFileIndex++);
    ofstream file(filename, std::ios::binary);
    if(!file.is_open()) {
        cout << "Could not open file " << filename << ". Aborting!" << endl;
        exit(1);
    }

    int numberOfPhaseSpaceCoordinates = 3*system->atoms().numberOfAtoms;
    float *phaseSpace = new float[numberOfPhaseSpaceCoordinates];
    vec3 systemSize = system->systemSize();

    int phaseSpaceCounter = 0;
    unsigned int numberOfAtoms = system->atoms().numberOfAtoms;

    for(int i=0;i<numberOfAtoms;i++) {
        phaseSpace[phaseSpaceCounter++] = system->atoms().x[i];
        phaseSpace[phaseSpaceCounter++] = system->atoms().y[i];
        phaseSpace[phaseSpaceCounter++] = system->atoms().z[i];
    }
    file.write((char*)&numberOfAtoms, sizeof(unsigned int));
    file.write((char*)(phaseSpace), numberOfPhaseSpaceCoordinates*sizeof(float));
    file.write((char*)(&systemSize[0]), 3*sizeof(MDDataType_t));

    file.close();
    delete phaseSpace;
}

void IO::savePositionsXYZ(System *system) {
    sprintf(filename, "state.xyz");
    ofstream file(filename, std::ios::binary);
    if(!file.is_open()) {
        cout << "Could not open file " << filename << ". Aborting!" << endl;
        exit(1);
    }

    file << system->atoms().numberOfAtoms << endl;
    file << "The is an optional comment line that can be empty." << endl;
    for(int n=0; n<system->atoms().numberOfAtoms; n++) {
        file << "Ar " << UnitConverter::lengthToAngstroms(system->atoms().x[n]) << " " << UnitConverter::lengthToAngstroms(system->atoms().y[n]) << " " << UnitConverter::lengthToAngstroms(system->atoms().z[n]) << endl;
    }
    file.close();
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

void IO::writeStatistics(float time, double kineticEnergy, double potentialEnergy, double pressure, double temperature) {
    if(!m_statisticsFile.is_open()) {
        m_statisticsFile.open("statistics.txt");
        m_statisticsFile << "Time [ps]              Kinetic energy [eV]    Potential energy [eV]   Total energy [eV]       Pressure [GPa]         Temperature [K]" << endl;
    }
    m_statisticsFile.precision(12);
    m_statisticsFile << std::scientific << UC::timeToSI(time)*1e12 << "     " << UC::energyToEv(kineticEnergy) << "     " << UC::energyToEv(potentialEnergy) << "     " << UC::energyToEv(kineticEnergy+potentialEnergy) << "     " << UC::pressureToSI(pressure)*1e-9 << "     " << UC::temperatureToSI(temperature) << endl;
}
