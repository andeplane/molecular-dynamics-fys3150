#include <io.h>
#include <system.h>
#include <atom.h>

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
    file << "The required line I have no idea why needs to be here." << endl;
    for(Atom *atom : system->atoms()) {
        file << "Ar " << atom->position.x << " " << atom->position.y << " " << atom->position.z << endl;
    }
}
