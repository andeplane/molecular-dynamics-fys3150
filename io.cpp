#include <io.h>
#include <system.h>
#include <atom.h>
#include <unitconverter.h>
using std::endl;

IO::IO()
{

}

void IO::open(char *filename) {
    file.open(filename);
}

void IO::close() {
    file.close();
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
void IO::saveState(System *system)
{
    file << system->atoms().size() << endl;
    file << "The is an optional comment line that can be empty." << endl;
    for(Atom *atom : system->atoms()) {
        file << "Ar " << UnitConverter::lengthToAngstroms(atom->position.x) << " " << UnitConverter::lengthToAngstroms(atom->position.y) << " " << UnitConverter::lengthToAngstroms(atom->position.z) << endl;
    }
}
