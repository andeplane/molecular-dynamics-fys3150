#ifndef CPELAPSEDTIMER_H
#define CPELAPSEDTIMER_H
#include <QElapsedTimer>
#include <QDebug>
class CPTimingObject {
private:
    QElapsedTimer m_timer;
    double m_timeElapsed;
public:
    CPTimingObject() : m_timeElapsed(0) { }

    void start() {
        m_timer.restart();
    }

    void stop() {
        m_timeElapsed += m_timer.elapsed() / double(1000);
        m_timer.restart();
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

    QElapsedTimer  m_timer;
    CPTimingObject m_calculateForces;
    CPTimingObject m_updateCellList;
    CPTimingObject m_updateNeighborList;
    CPTimingObject m_move;
    CPTimingObject m_halfKick;
    CPTimingObject m_periodicBoundaryConditions;
    CPTimingObject m_sampling;
    CPTimingObject m_disk;

    static CPTimingObject &calculateForces() { return CPElapsedTimer::getInstance().m_calculateForces; }
    static CPTimingObject &updateCellList() { return CPElapsedTimer::getInstance().m_updateCellList; }
    static CPTimingObject &updateNeighborList() { return CPElapsedTimer::getInstance().m_updateNeighborList; }
    static CPTimingObject &move() { return CPElapsedTimer::getInstance().m_move; }
    static CPTimingObject &halfKick() { return CPElapsedTimer::getInstance().m_halfKick; }
    static CPTimingObject &periodicBoundaryConditions() { return CPElapsedTimer::getInstance().m_periodicBoundaryConditions; }
    static CPTimingObject &sampling() { return CPElapsedTimer::getInstance().m_sampling; }
    static CPTimingObject &disk() { return CPElapsedTimer::getInstance().m_disk; }
    static double totalTime() { return CPElapsedTimer::getInstance().m_timer.elapsed() / double(1000); }
};

#endif // CPELAPSEDTIMER_H
