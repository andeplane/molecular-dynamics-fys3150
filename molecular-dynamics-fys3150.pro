TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
#CONFIG -= qt
INCLUDEPATH +=
# QT += quick widgets
CONFIG += c++11
#QMAKE_CXX = g++-4.8

QMAKE_LINK = icpc
QMAKE_CC = icc
QMAKE_CXX = icpc
QMAKE_CXXFLAGS += -O3
QMAKE_CFLAGS += -O3
#QMAKE_CXXFLAGS += -xCORE-AVX-I -O3 -ipo -g -falign-functions=16
#QMAKE_CFLAGS += -xCORE-AVX-I -O3 -ipo -g -falign-functions=16

#QMAKE_CXXFLAGS += -prof-gen
#QMAKE_CFLAGS += -prof-gen

#QMAKE_CXXFLAGS += -prof-use
#QMAKE_CFLAGS += -prof-use

DEFINES += MD_SIMD
#DEFINES += MD_DEBUG

SOURCES += main.cpp \
    system.cpp \
    integrators/integrator.cpp \
    integrators/velocityverlet.cpp \
    math/vec3.cpp \
    math/random.cpp \
    io.cpp \
    potentials/potential.cpp \
    potentials/lennardjones.cpp \
    statisticssampler.cpp \
    unitconverter.cpp \
    celllist.cpp \
    neighborlist.cpp \
    cpelapsedtimer.cpp \
    modifiers/berendsenthermostat.cpp \
    math/hilbert.cpp \
    math/morton.cpp \
    atoms.cpp

HEADERS += \
    system.h \
    integrators/integrator.h \
    integrators/velocityverlet.h \
    math/vec3.h \
    math/random.h \
    io.h \
    potentials/potential.h \
    potentials/lennardjones.h \
    statisticssampler.h \
    unitconverter.h \
    celllist.h \
    neighborlist.h \
    cpelapsedtimer.h \
    modifiers/berendsenthermostat.h \
    math/hilbert.h \
    math/morton.h \
    atoms.h \
    config.h

