#pragma once
#include <fstream>
class System;
using std::ofstream;

class IO
{
private:
    ofstream file;
public:
    IO();
    ~IO();

    void saveState(System *system);
    void open(char *filename);
    void close();

};
