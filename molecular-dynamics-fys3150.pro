TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    atom.cpp \
    system.cpp \
    integrators/integrator.cpp \
    integrators/velocityverlet.cpp \
    math/vec3.cpp \
    math/random.cpp

HEADERS += \
    atom.h \
    system.h \
    integrators/integrator.h \
    integrators/velocityverlet.h \
    math/vec3.h \
    math/random.h

