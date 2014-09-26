#include <math/vec3.h>
#include <cmath>
namespace CompPhys {
    vec3::vec3() :
        x(0),
        y(0),
        z(0)
    {

    }

    vec3::vec3(double x_, double y_, double z_) :
        x(x_),
        y(y_),
        z(z_)
    {

    }

    bool vec3::operator==(vec3 &rhs) {
        return(x == rhs.x && y == rhs.y && z == rhs.z);
    }

    vec3 vec3::operator+(vec3 &rhs) {
        return vec3( x + rhs.x,
                     y + rhs.y,
                     z + rhs.z);
    }

    vec3 vec3::operator-(vec3 &rhs) {
        return vec3( x - rhs.x,
                     y - rhs.y,
                     z - rhs.z);
    }

    vec3 vec3::operator*(vec3 &rhs) {
        return vec3( x * rhs.x,
                     y * rhs.y,
                     z * rhs.z);
    }

    vec3 vec3::operator/(vec3 &rhs) {
        return vec3( x / rhs.x,
                     y / rhs.y,
                     z / rhs.z);
    }

    vec3 vec3::operator+(double scalar) {
        return vec3(x + scalar,
                    y + scalar,
                    z + scalar);
    }

    vec3 vec3::operator-(double scalar) {
        return vec3(x - scalar,
                    y - scalar,
                    z - scalar);
    }

    vec3 vec3::operator*(double scalar) {
        return vec3(x * scalar,
                    y * scalar,
                    z * scalar);
    }

    vec3 vec3::operator/(double scalar) {
        return vec3(x / scalar,
                    y / scalar,
                    z / scalar);
    }

    void vec3::addAndMultiply(vec3 &rhs, double scalar) {
        x += rhs.x*scalar;
        y += rhs.y*scalar;
        z += rhs.z*scalar;
    }

    void vec3::add(vec3 &rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
    }

    double vec3::dot(vec3 &rhs) {
        return (x * rhs.x +
                y * rhs.y +
                z * rhs.z);
    }

    vec3 vec3::cross(vec3 &rhs) {
        return vec3( y * rhs.z - z * rhs.y,
                     z * rhs.x - x * rhs.z,
                     x * rhs.y - y * rhs.x);
    }

    double vec3::length() {
        return sqrt(lengthSquared());
    }

    double vec3::lengthSquared() {
        return x*x + y*y + z*z;
    }

    void vec3::normalize() {
        double myLength = length();
        if(myLength > 0) { // Don't divide by zero...
            x /= myLength;
            y /= myLength;
            z /= myLength;
        }
    }

    void vec3::setToZero()
    {
        x = 0;
        y = 0;
        z = 0;
    }
}
