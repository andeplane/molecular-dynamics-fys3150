TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG += c++11
CONFIG -= qt
unix:!macx {
    QMAKE_CXXFLAGS += -I/usr/include/openmpi-x86_64 -pthread -m64
    LIBS += -pthread -m64 -Wl,-rpath -Wl,/usr/lib64/openmpi/lib -Wl,--enable-new-dtags -L/usr/lib64/openmpi/lib -lmpi_cxx -lmpi

}

macx: {
    INCLUDEPATH += /usr/local/Cellar/open-mpi/2.0.1/include
    LIBS += -L/usr/local/opt/libevent/lib -L/usr/local/Cellar/open-mpi/2.0.1/lib -lmpi
}


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
    math/morton.cpp

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
    math/morton.h

