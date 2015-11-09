#include "vec3.h"
#include "random.h"
#include <cmath>
using std::cout; using std::endl;

vec3::vec3()
{
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
}

vec3::vec3(const vec3 &copy)
{
    components[0] = copy.x();
    components[1] = copy.y();
    components[2] = copy.z();
}

vec3::vec3(double x, double y, double z)
{
    components[0] = x;
    components[1] = y;
    components[2] = z;
}

void vec3::print()
{
    // Will print matlab syntax vector. Output will be like: [2.09, 5.3, 9.1];
    cout << "[" << components[0] << ", " << components[1] << ", " << components[2] << "]" << endl;
}

void vec3::print(string name)
{
    // Will print matlab syntax vector with a name. Output will be like: A = [2.09, 5.3, 9.1];
    cout << name << " = ";
    print();
}

vec3 vec3::cross(vec3 aVector)
{
    return vec3(y()*aVector.z()-z()*aVector.y(), z()*aVector.x()-x()*aVector.z(), x()*aVector.y()-y()*aVector.x());
}

double vec3::lengthSquared()
{
    // Returns the square of the length (or norm) of the vector
    return components[0]*components[0]+components[1]*components[1]+components[2]*components[2];
}

double vec3::length()
{
    // Returns the length (or norm) of the vector
    return sqrt(lengthSquared());
}

void vec3::set(double x, double y, double z)
{
    setX(x);
    setY(y);
    setZ(z);
}

void vec3::randomGaussian(double mean, double standardDeviation)
{
    components[0] = Random::nextGaussian(mean, standardDeviation);
    components[1] = Random::nextGaussian(mean, standardDeviation);
    components[2] = Random::nextGaussian(mean, standardDeviation);
}

vec3 &vec3::operator+=(double rhs)
{
    components[0] += rhs;
    components[1] += rhs;
    components[2] += rhs;
    return *this;
}

vec3 &vec3::operator+=(vec3 rhs)
{
    components[0] += rhs[0];
    components[1] += rhs[1];
    components[2] += rhs[2];
    return *this;
}

vec3 &vec3::operator*=(double rhs)
{
    components[0] *= rhs;
    components[1] *= rhs;
    components[2] *= rhs;
    return *this;
}

vec3 &vec3::operator*=(vec3 rhs)
{
    components[0] *= rhs[0];
    components[1] *= rhs[1];
    components[2] *= rhs[2];
    return *this;
}

vec3 &vec3::operator-=(double rhs)
{
    components[0] -= rhs;
    components[1] -= rhs;
    components[2] -= rhs;
    return *this;
}

vec3 &vec3::operator-=(vec3 rhs)
{
    components[0] -= rhs[0];
    components[1] -= rhs[1];
    components[2] -= rhs[2];
    return *this;
}

vec3 &vec3::operator/=(double rhs)
{
    components[0] /= rhs;
    components[1] /= rhs;
    components[2] /= rhs;
    return *this;
}

vec3 &vec3::operator/=(vec3 rhs)
{
    components[0] /= rhs[0];
    components[1] /= rhs[1];
    components[2] /= rhs[2];
    return *this;
}

std::ostream &operator<<(std::ostream &os, const vec3 &myVector) // Allows cout << myVector << endl;
{
    os << "[" << myVector.x() << ", " << myVector.y() << ", " << myVector.z() << "];";
    return os;
}
