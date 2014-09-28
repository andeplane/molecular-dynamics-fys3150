#pragma once
#include <iostream>

namespace CompPhys {
    class vec3
    {
    private:
        friend std::ostream& operator<<(std::ostream&stream, vec3 &vec);
    public:
        vec3(); // Create a zero vector
        vec3(double x, double y, double z);
        bool operator==(vec3 &rhs);
        vec3 operator+(vec3 &rhs);
        vec3 operator-(vec3 &rhs);
        vec3 operator*(vec3 &rhs);
        vec3 operator/(vec3 &rhs);
        vec3 operator+(double scalar);
        vec3 operator-(double scalar);
        vec3 operator*(double scalar);
        vec3 operator/(double scalar);
        void add(vec3 &rhs);
        void addAndMultiply(vec3 &rhs, double scalar);
        vec3 cross(vec3 &rhs);
        double dot(vec3 &rhs);
        double length();
        double lengthSquared();
        void normalize();
        void setToZero();

        double x;
        double y;
        double z;
    };
}
