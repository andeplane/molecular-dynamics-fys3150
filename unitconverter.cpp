#include <iostream>
#include <cmath>
#include "unitconverter.h"
#include "math/vec3.h"

using namespace std;

float UnitConverter::m0 = 0;
float UnitConverter::q0 = 0;
float UnitConverter::hbar0 = 0;
float UnitConverter::electricConstant0 = 0;
float UnitConverter::a0 = 0;
float UnitConverter::a0_angstrom = 0;
float UnitConverter::E0 = 0;
float UnitConverter::E0ev = 0;
float UnitConverter::kb = 0;
float UnitConverter::t0 = 0;
float UnitConverter::F0 = 0;
float UnitConverter::T0 = 0;
float UnitConverter::P0 = 0;
float UnitConverter::v0 = 0;
float UnitConverter::visc0 = 0;
float UnitConverter::diff0 = 0;
std::string UnitConverter::currentUnits = "No units chosen.";

bool UnitConverter::initialized = false;

void UnitConverter::initializeMDUnits() {
    UnitConverter::initialized = true;
    UnitConverter::currentUnits = "MD units";
    // Molecular Dynamics units

    // Fundamental units
    float m0 = 1.66053892e-27;         // SI [kg]
    float L0 = 1e-10;                  // SI [m]
    float kb = 1.3806488e-23;          // SI [J/K]
    float E0eV = 1.0318e-2;            // eV
    float E0 = 1.60217657e-19*E0eV;    // SI [J]

    UnitConverter::m0 = m0;
    UnitConverter::kb = kb;

    // Derived units
    UnitConverter::a0 = L0;
    UnitConverter::E0 = E0;
    UnitConverter::t0 = L0*sqrt(m0/E0);
    UnitConverter::v0 = L0/t0;
    UnitConverter::F0 = E0/a0;
    UnitConverter::T0 = E0/kb;
    UnitConverter::P0 = E0/(a0*a0*a0);
    UnitConverter::visc0 = P0*t0;
    UnitConverter::diff0 = a0*a0/t0;
    UnitConverter::E0ev = E0eV;
    UnitConverter::hbar0 = INFINITY; // Not really used, give INF as a warning-ish?
    UnitConverter::q0 = INFINITY;
    UnitConverter::electricConstant0 = INFINITY;
}

void UnitConverter::initializeAtomicUnits() {
    UnitConverter::initialized = true;
    UnitConverter::currentUnits = "Atomic units";
    // Atomic units
    // [1] http://en.wikipedia.org/wiki/Atomic_units

    // Fundamental units
    float m0 = 9.10938291e-31;  // SI [kg]
    float q0 = 1.602176565e-19; // SI [C]
    float hbar0 = 1.054571726e-34; // SI [Js]
    float electricConstant0 = 8.9875517873681e9; // SI [kgm^3/(s^-2 C^-2)]
    float kb = 1.3806488e-23; // SI [J/K]

    UnitConverter::m0 = m0;
    UnitConverter::q0 = q0;
    UnitConverter::hbar0 = hbar0;
    UnitConverter::electricConstant0 = electricConstant0;
    UnitConverter::kb = kb;

    // Derived units
    UnitConverter::a0 = hbar0*hbar0/(electricConstant0*m0*q0*q0);
    UnitConverter::E0 = m0*q0*q0*q0*q0*electricConstant0*electricConstant0/(hbar0*hbar0);
    UnitConverter::t0 = hbar0/E0;
    UnitConverter::v0 = a0*E0/hbar0;
    UnitConverter::F0 = E0/a0;
    UnitConverter::T0 = E0/kb;
    UnitConverter::P0 = E0/(a0*a0*a0);
    UnitConverter::visc0 = P0*t0;
    UnitConverter::diff0 = a0*a0/t0;

    UnitConverter::E0ev = 1.0/UnitConverter::energyFromSI(1.60217657e-19);
}

void UnitConverter::makeSureInitialized() {
    if(!UnitConverter::initialized) UnitConverter::initialize(MDUnits);
}

void UnitConverter::initialize(Units type) {
    if(type == AtomicUnits) UnitConverter::initializeAtomicUnits();
    else if(type == MDUnits) UnitConverter::initializeMDUnits();
}

float UnitConverter::pressureToSI(float P) {UnitConverter::makeSureInitialized(); return UnitConverter::P0*P; }
float UnitConverter::pressureFromSI(float P) {UnitConverter::makeSureInitialized(); return P/UnitConverter::P0; }

float UnitConverter::temperatureToSI(float T) {UnitConverter::makeSureInitialized(); return UnitConverter::T0*T; }
float UnitConverter::temperatureFromSI(float T) {UnitConverter::makeSureInitialized(); return T/UnitConverter::T0; }

float UnitConverter::massToSI(float m) {UnitConverter::makeSureInitialized(); return UnitConverter::m0*m; }
float UnitConverter::massFromSI(float m) {UnitConverter::makeSureInitialized(); return m/UnitConverter::m0; }

float UnitConverter::lengthToSI(float L) {UnitConverter::makeSureInitialized(); return UnitConverter::a0*L; }
float UnitConverter::lengthFromSI(float L) {UnitConverter::makeSureInitialized(); return L/UnitConverter::a0; }

float UnitConverter::lengthToAngstroms(float L) {UnitConverter::makeSureInitialized(); return UnitConverter::a0*L*1e10; }
float UnitConverter::lengthFromAngstroms(float L) {UnitConverter::makeSureInitialized(); return L/(UnitConverter::a0*1e10); }

vec3 UnitConverter::lengthToSI(vec3 position)
{
    return vec3(UnitConverter::lengthToSI(position.x()), UnitConverter::lengthToSI(position.y()), UnitConverter::lengthToSI(position.z()));
}

vec3 UnitConverter::lengthFromSI(vec3 position)
{
    return vec3(UnitConverter::lengthFromSI(position.x()), UnitConverter::lengthFromSI(position.y()), UnitConverter::lengthFromSI(position.z()));
}

vec3 UnitConverter::lengthToAngstroms(vec3 position)
{
    return vec3(UnitConverter::lengthToAngstroms(position.x()), UnitConverter::lengthToAngstroms(position.y()), UnitConverter::lengthToAngstroms(position.z()));
}

vec3 UnitConverter::lengthFromAngstroms(vec3 position)
{
    return vec3(UnitConverter::lengthFromAngstroms(position.x()), UnitConverter::lengthFromAngstroms(position.y()), UnitConverter::lengthFromAngstroms(position.z()));
}

vec3 UnitConverter::velocityToSI(vec3 velocity)
{
    return vec3(UnitConverter::velocityToSI(velocity.x()), UnitConverter::velocityToSI(velocity.y()), UnitConverter::velocityToSI(velocity.z()));
}

vec3 UnitConverter::velocityFromSI(vec3 velocity)
{
    return vec3(UnitConverter::velocityFromSI(velocity.x()), UnitConverter::velocityFromSI(velocity.y()), UnitConverter::velocityFromSI(velocity.z()));
}

float UnitConverter::forceToSI(float F) {UnitConverter::makeSureInitialized(); return UnitConverter::F0*F; }
float UnitConverter::forceFromSI(float F) {UnitConverter::makeSureInitialized(); return F/UnitConverter::F0; }

float UnitConverter::energyToSI(float E) {UnitConverter::makeSureInitialized(); return UnitConverter::E0*E; }
float UnitConverter::energyFromSI(float E) {UnitConverter::makeSureInitialized(); return E/UnitConverter::E0; }

float UnitConverter::energyToEv(float E) {UnitConverter::makeSureInitialized(); return UnitConverter::E0ev*E; }
float UnitConverter::energyFromEv(float E) {UnitConverter::makeSureInitialized(); return E/UnitConverter::E0ev; }

float UnitConverter::degreesToRadians(float v) {UnitConverter::makeSureInitialized(); return M_PI/180*v; }
float UnitConverter::radiansToDegrees(float v) {UnitConverter::makeSureInitialized(); return 180/M_PI*v; }

float UnitConverter::timeToSI(float t) {UnitConverter::makeSureInitialized(); return UnitConverter::t0*t; }
float UnitConverter::timeFromSI(float t) {UnitConverter::makeSureInitialized(); return t/UnitConverter::t0; }

float UnitConverter::velocityToSI(float v) {UnitConverter::makeSureInitialized(); return v*UnitConverter::v0; }
float UnitConverter::velocityFromSI(float v) {UnitConverter::makeSureInitialized(); return v/UnitConverter::v0; }

float UnitConverter::diffusionToSI(float d) {UnitConverter::makeSureInitialized(); return d*UnitConverter::diff0; }
float UnitConverter::diffusionFromSI(float d) {UnitConverter::makeSureInitialized(); return d/UnitConverter::diff0; }
