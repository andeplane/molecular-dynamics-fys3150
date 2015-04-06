#include <cmath>
#include <atom.h>
#include "math/random.h"
#include <iostream>
#include <cassert>

using namespace std;

int Atom::totalNumberOfAtoms = 0;


unsigned int Atom::cellIndex() const
{
    return m_cellIndex;
}

void Atom::setCellIndex(unsigned int cellIndex)
{
    m_cellIndex = cellIndex;
}
Atom::Atom() :
    m_mass(NAN),
    m_index(Atom::totalNumberOfAtoms++),
    m_cellIndex(0)
{

}

Atom::~Atom()
{

}

void Atom::resetForce()
{
    force.setToZero();
}

void Atom::resetVelocityMaxwellian(float temperature)
{
    // Resetting the velocity according to a Maxwell-Boltzmann distribution (see http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution )
    float boltzmannConstant = 1.0; // In atomic units, the boltzmann constant equals 1
    float standardDeviation = sqrt(boltzmannConstant*temperature/m_mass);
    velocity.randomGaussian(0, standardDeviation);
}

void Atom::resetNeighbors() {
    // m_numNeighbors = 0;
    // m_extraNeighbors.clear();
}

void Atom::addNeighbor(Atom *atom)
{
//    if(m_numNeighbors>=50) {
//        m_extraNeighbors.push_back(atom);
//    } else {
//        m_neighbors[m_numNeighbors++] = atom;
//    }
    // assert(m_numNeighbors < 150 && "Too many neighbors...");

    // m_neighbors[m_numNeighbors++] = atom;
}
