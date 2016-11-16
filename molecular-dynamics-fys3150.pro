TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX += -g

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    velocityverlet.cpp \
    math/vec3.cpp \
    io.cpp \
    lennardjones.cpp \
    statisticssampler.cpp \
    unitconverter.cpp

HEADERS += \
    atom.h \
    system.h \
    velocityverlet.h \
    math/vec3.h \
    math/random.h \
    io.h \
    lennardjones.h \
    statisticssampler.h \
    unitconverter.h

