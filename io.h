#pragma once
#include <fstream>
class System;
using std::ofstream;

class IO
{
private:
    ofstream m_movieFile;
    ofstream m_statisticsFile;
    ofstream m_performanceFile;
public:
    IO();
    ~IO();

    void saveState(System *system);
    void close();

    void writeStatistics(float time, float kineticEnergy, float potentialEnergy, float pressure, float temperature);
    void writePerformance(unsigned int timestep);
};
