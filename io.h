#pragma once
#include <fstream>
class System;
using std::ofstream;

class IO
{
public:
    IO();
    ofstream file;
    void saveState(System *system);
    void open(char *filename);
    void close();
};
