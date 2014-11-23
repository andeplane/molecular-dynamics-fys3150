TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
#CONFIG -= qt
INCLUDEPATH +=
# QT += quick widgets
CONFIG += c++11
QMAKE_CXX = g++-4.9

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    integrators/integrator.cpp \
    integrators/velocityverlet.cpp \
    math/vec3.cpp \
    math/random.cpp \
    io.cpp \
    potentials/potential.cpp \
    potentials/lennardjones.cpp \
    statisticssampler.cpp \
    integrators/eulercromer.cpp \
    unitconverter.cpp \
    celllist.cpp \
    neighborlist.cpp \
    cpelapsedtimer.cpp \
    modifiers/berendsenthermostat.cpp \
    math/hilbert.cpp \
    math/morton.cpp \
    atoms.cpp

HEADERS += \
    atom.h \
    system.h \
    integrators/integrator.h \
    integrators/velocityverlet.h \
    math/vec3.h \
    math/random.h \
    io.h \
    potentials/potential.h \
    potentials/lennardjones.h \
    statisticssampler.h \
    integrators/eulercromer.h \
    unitconverter.h \
    celllist.h \
    neighborlist.h \
    cpelapsedtimer.h \
    modifiers/berendsenthermostat.h \
    math/hilbert.h \
    math/morton.h \
    atoms.h

