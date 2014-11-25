#include "cpelapsedtimer.h"

CPElapsedTimer::CPElapsedTimer()
{
    m_startedAt = clock();
    m_lastPing = m_startedAt;
}

double CPElapsedTimer::ping()
{
    CPElapsedTimer &timer = CPElapsedTimer::getInstance();

    double deltaT = CPElapsedTimer::totalTime() - timer.m_lastPing;
    timer.m_lastPing = CPElapsedTimer::totalTime();
    return deltaT;
}
