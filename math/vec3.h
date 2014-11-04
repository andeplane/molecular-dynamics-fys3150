#pragma once
#include <iostream>

class vec3
{
private:
    float m_vec[3];
public:
    vec3(); // Create a zero vector
    vec3(float x, float y, float z);
    bool operator==(vec3 rhs);
    vec3 operator+(vec3 rhs);
    vec3 &operator+=(vec3 rhs);
    vec3 operator-(vec3 rhs);
    vec3 &operator-=(vec3 rhs);
    vec3 operator*(vec3 rhs);
    vec3 &operator*=(vec3 rhs);
    vec3 operator/(vec3 rhs);
    vec3 &operator/=(vec3 rhs);
    vec3 operator+(float scalar);
    vec3 &operator+=(float scalar);
    vec3 operator-(float scalar);
    vec3 &operator-=(float scalar);
    vec3 operator*(float scalar);
    vec3 &operator*=(float scalar);
    vec3 operator/(float scalar);
    vec3 &operator/=(float scalar);
    vec3 operator-();
    void add(vec3 &rhs) {
        m_vec[0] += rhs.x();
        m_vec[1] += rhs.y();
        m_vec[2] += rhs.z();
    }
    void addAndMultiply(vec3 &rhs, float scalar) {
        m_vec[0] += rhs.x()*scalar;
        m_vec[1] += rhs.y()*scalar;
        m_vec[2] += rhs.z()*scalar;
    }
    vec3 cross(vec3 &rhs);
    float dot(vec3 &rhs);
    float length();
    float lengthSquared();
    void normalize();
    void setToZero();
    void randomGaussian(float mean, float standardDeviation);
    void randomUniform(float min, float max);
    void set(float x, float y, float z);
    inline float x() const { return m_vec[0]; }
    inline float y() const { return m_vec[1]; }
    inline float z() const { return m_vec[2]; }
    inline float &operator[](int index) { return m_vec[index]; }
    inline float operator[](int index) const { return m_vec[index]; }
private:
    friend std::ostream& operator<<(std::ostream&stream, vec3 vec);
};
