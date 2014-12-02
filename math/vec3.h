#pragma once
#include <iostream>
#include "config.h"

class vec3
{
private:
    MDDataType_t m_vec[3];
public:
    vec3(); // Create a zero vector
    vec3(MDDataType_t x, MDDataType_t y, MDDataType_t z);
    bool operator==(vec3 rhs);
    vec3 operator+(vec3 rhs);
    vec3 &operator+=(vec3 rhs);
    vec3 operator-(vec3 rhs);
    vec3 &operator-=(vec3 rhs);
    vec3 operator*(vec3 rhs);
    vec3 &operator*=(vec3 rhs);
    vec3 operator/(vec3 rhs);
    vec3 &operator/=(vec3 rhs);
    vec3 operator+(MDDataType_t scalar);
    vec3 &operator+=(MDDataType_t scalar);
    vec3 operator-(MDDataType_t scalar);
    vec3 &operator-=(MDDataType_t scalar);
    vec3 operator*(MDDataType_t scalar);
    vec3 &operator*=(MDDataType_t scalar);
    vec3 operator/(MDDataType_t scalar);
    vec3 &operator/=(MDDataType_t scalar);
    inline vec3 operator-() { return vec3(-m_vec[0], -m_vec[1], -m_vec[2]); }
    void add(vec3 &rhs) {
        m_vec[0] += rhs.x();
        m_vec[1] += rhs.y();
        m_vec[2] += rhs.z();
    }
    void addAndMultiply(vec3 &rhs, MDDataType_t scalar) {
        m_vec[0] += rhs.x()*scalar;
        m_vec[1] += rhs.y()*scalar;
        m_vec[2] += rhs.z()*scalar;
    }
    vec3 cross(vec3 &rhs);
    MDDataType_t dot(vec3 &rhs);
    MDDataType_t length();
    void normalize();
    void setToZero();
    void randomGaussian(MDDataType_t mean, MDDataType_t standardDeviation);
    void randomUniform(MDDataType_t min, MDDataType_t max);
    void set(MDDataType_t x, MDDataType_t y, MDDataType_t z);
    inline MDDataType_t x() const { return m_vec[0]; }
    inline MDDataType_t y() const { return m_vec[1]; }
    inline MDDataType_t z() const { return m_vec[2]; }
    inline MDDataType_t &operator[](int index) { return m_vec[index]; }
    inline MDDataType_t operator[](int index) const { return m_vec[index]; }
    inline MDDataType_t lengthSquared() { return m_vec[0]*m_vec[0] + m_vec[1]*m_vec[1] + m_vec[2]*m_vec[2]; }
    inline void subtract(const vec3 &v1, const vec3 &v2) { m_vec[0] = v1[0] - v2[0]; m_vec[1] = v1[1] - v2[1]; m_vec[2] = v1[2] - v2[2]; }
private:
    friend std::ostream& operator<<(std::ostream&stream, vec3 vec);
};
