#include <io.h>
#include <system.h>
#include <atom.h>
#include <unitconverter.h>
#include <cstdlib>
using std::endl; using std::cout;

IO::IO()
{

}

IO::~IO() {
    close();
}

void IO::open(char *filename) {
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
void IO::saveState(System *system)
{
    file << system->atoms().size() << endl;
    file << "The is an optional comment line that can be empty." << endl;
    for(int n=0; n<system->atoms().size(); n++) {
        Atom *atom = system->atoms()[n];
        file << "Ar " << UnitConverter::lengthToAngstroms(atom->position.x()) << " " << UnitConverter::lengthToAngstroms(atom->position.y()) << " " << UnitConverter::lengthToAngstroms(atom->position.z()) << endl;
    }
}
