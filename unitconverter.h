#pragma once
#include <string>
#include <vector>
#include <math/vec3.h>

using std::vector;
using CompPhys::vec3;

enum Units {SIUnits = 0, AtomicUnits, MDUnits};

class UnitConverter
{
public:
    static double m0;
    static double q0;
    static double hbar0;
    static double electricConstant0;
    static double a0;
    static double a0_angstrom;
    static double E0;
    static double E0ev;
    static double kb;
    static double t0;
    static double F0;
    static double T0;
    static double P0;
    static double v0;
    static double visc0;
    static double diff0;
    static std::string currentUnits;

    static void makeSureInitialized();
    static void initialize(Units type);
    static bool initialized;

    static double pressureToSI(double P);
    static double pressureFromSI(double P);

    static double temperatureToSI(double T);
    static double temperatureFromSI(double T);

    static double massToSI(double m);
    static double massFromSI(double m);

    static double lengthToSI(double L);
    static double lengthFromSI(double L);
    static vec3 lengthToSI(vec3 position);
    static vec3 lengthFromSI(vec3 position);

    static vec3 velocityToSI(vec3 position);
    static vec3 velocityFromSI(vec3 position);

    static vec3 lengthToAngstroms(vec3 position);
    static vec3 lengthFromAngstroms(vec3 position);

    static double lengthToAngstroms(double L);
    static double lengthFromAngstroms(double L);

    static double forceToSI(double F);
    static double forceFromSI(double F);

    static double energyToSI(double E);
    static double energyFromSI(double E);

    static double timeToSI(double t);
    static double timeFromSI(double t);

    static double velocityToSI(double v);
    static double velocityFromSI(double v);

    static double diffusionToSI(double d);
    static double diffusionFromSI(double d);

    static double energyToEv(double E);
    static double energyFromEv(double E);

    static double degreesToRadians(double v);
    static double radiansToDegrees(double v);

    static void initializeAtomicUnits();
    static void initializeMDUnits();
};
