#include <fstream>
#include <io.h>
#include <system.h>
#include <atom.h>

using std::ofstream; using std::endl;

IO::IO()
{

}

void IO::saveStateToFile(System *system, char *filename)
{
    ofstream file(filename);

    for(Atom *atom : system->atoms()) {
        file << "Ar " << atom->position.x << " " << atom->position.y << " " << atom->position.z << endl;
    }

    file.close();
}
