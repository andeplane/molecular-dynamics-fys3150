#ifndef CPELAPSEDTIMER_H
#define CPELAPSEDTIMER_H
#include <time.h>

class CPTimingObject {
private:
    double m_timeElapsed;
    clock_t m_startedAt;
public:
    CPTimingObject() : m_timeElapsed(0), m_startedAt(0) { }

    void start() {
         m_startedAt = clock();
    }

    double stop() {
        double t = double(clock() - m_startedAt)/CLOCKS_PER_SEC;
        m_timeElapsed += t;
        return t;
    }

    double elapsedTime() { return m_timeElapsed; }
};

class CPElapsedTimer
{
public:
    CPElapsedTimer();

    static CPElapsedTimer& getInstance()
    {
        static CPElapsedTimer instance; // Guaranteed to be destroyed.
                                 // Instantiated on first use.
        return instance;
    }

    clock_t        m_startedAt;
    CPTimingObject m_calculateForces;
    CPTimingObject m_updateCellList;
    CPTimingObject m_updateNeighborList;
    CPTimingObject m_updateNeighborCopies;
    CPTimingObject m_move;
    CPTimingObject m_halfKick;
    CPTimingObject m_periodicBoundaryConditions;
    CPTimingObject m_sampling;
    CPTimingObject m_disk;
    CPTimingObject m_timeEvolution;
    CPTimingObject m_thermostat;

    static CPTimingObject &calculateForces() { return CPElapsedTimer::getInstance().m_calculateForces; }
    static CPTimingObject &updateCellList() { return CPElapsedTimer::getInstance().m_updateCellList; }
    static CPTimingObject &updateNeighborList() { return CPElapsedTimer::getInstance().m_updateNeighborList; }
    static CPTimingObject &updateNeighborCopies() { return CPElapsedTimer::getInstance().m_updateNeighborCopies; }
    static CPTimingObject &move() { return CPElapsedTimer::getInstance().m_move; }
    static CPTimingObject &halfKick() { return CPElapsedTimer::getInstance().m_halfKick; }
    static CPTimingObject &periodicBoundaryConditions() { return CPElapsedTimer::getInstance().m_periodicBoundaryConditions; }
    static CPTimingObject &sampling() { return CPElapsedTimer::getInstance().m_sampling; }
    static CPTimingObject &disk() { return CPElapsedTimer::getInstance().m_disk; }
    static CPTimingObject &timeEvolution() { return CPElapsedTimer::getInstance().m_timeEvolution; }
    static CPTimingObject &thermostat() { return CPElapsedTimer::getInstance().m_thermostat; }

    static double totalTime() { return double(clock() - CPElapsedTimer::getInstance().m_startedAt)/ CLOCKS_PER_SEC; }
};

#endif // CPELAPSEDTIMER_H
