#include "io.h"
#include "system.h"
#include "atom.h"
#include "unitconverter.h"
#include <cstdlib>
using std::endl; using std::cout;


IO::IO(const char *filename)
{
    open(filename);
}

IO::~IO() {
    close();
}

void IO::open(const char *filename) {
    if(file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
        exit(1);
    }

    file.open(filename);
}

void IO::close() {
    if(file.is_open()) {
        file.close();
    }
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
// It can easily be opened in Ovito. Note that you can also output more properties than just the position. You can print the
// velocities per particle (or kinetic energy etc), and color the atoms in Ovito based on these properties.

void IO::saveState(System &system)
{
    if(file.is_open()) {
        file << system.atoms().size() << endl;
        file << "The is an optional comment line that can be empty. The reason we use H is so particles get smaller in Ovito" << endl;
        for(Atom *atom : system.atoms()) {
            file << "H " <<
                    UnitConverter::lengthToAngstroms(atom->position.x()) << " " <<
                    UnitConverter::lengthToAngstroms(atom->position.y()) << " " <<
                    UnitConverter::lengthToAngstroms(atom->position.z()) << "\n";
        }
    }
}
