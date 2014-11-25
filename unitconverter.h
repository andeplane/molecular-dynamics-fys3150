#pragma once
#include <string>
#include <vector>
#include <math/vec3.h>

using std::vector;

enum Units {SIUnits = 0, AtomicUnits, MDUnits};

class UnitConverter
{
public:
    static float m0;
    static float q0;
    static float hbar0;
    static float electricConstant0;
    static float a0;
    static float a0_angstrom;
    static float E0;
    static float E0ev;
    static float kb;
    static float t0;
    static float F0;
    static float T0;
    static float P0;
    static float v0;
    static float visc0;
    static float diff0;
    static std::string currentUnits;

    static void makeSureInitialized();
    static void initialize(Units type);
    static bool initialized;

    static float pressureToSI(float P);
    static float pressureFromSI(float P);

    static float temperatureToSI(float T);
    static float temperatureFromSI(float T);

    static float massToSI(float m);
    static float massFromSI(float m);

    static float lengthToSI(float L);
    static float lengthFromSI(float L);
    static vec3 lengthToSI(vec3 position);
    static vec3 lengthFromSI(vec3 position);

    static vec3 velocityToSI(vec3 position);
    static vec3 velocityFromSI(vec3 position);

    static vec3 lengthToAngstroms(vec3 position);
    static vec3 lengthFromAngstroms(vec3 position);

    static float lengthToAngstroms(float L);
    static float lengthFromAngstroms(float L);

    static float forceToSI(float F);
    static float forceFromSI(float F);

    static float energyToSI(float E);
    static float energyFromSI(float E);

    static float timeToSI(float t);
    static float timeFromSI(float t);

    static float velocityToSI(float v);
    static float velocityFromSI(float v);

    static float diffusionToSI(float d);
    static float diffusionFromSI(float d);

    static float energyToEv(float E);
    static float energyFromEv(float E);

    static float degreesToRadians(float v);
    static float radiansToDegrees(float v);

    static void initializeAtomicUnits();
    static void initializeMDUnits();
};

typedef UnitConverter UC;
