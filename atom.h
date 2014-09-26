#pragma once
#include <math/vec3.h>
using fys3150::vec3;

class Atom
{
public:
    vec3 position;
    vec3 velocity;
    vec3 force;

    Atom();
    void resetForce();
};
