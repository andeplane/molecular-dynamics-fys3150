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
    inline vec3 operator-() { return vec3(-m_vec[0], -m_vec[1], -m_vec[2]); }
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
    void normalize();
    void setToZero();
    void randomGaussian(float mean, float standardDeviation);
    void randomUniform(float min, float max);
    void set(float x, float y, float z);
    inline float x() const { return m_vec[0]; }
    inline float y() const { return m_vec[1]; }
    inline float z() const { return m_vec[2]; }
    inline void addX(float val) { m_vec[0] += val; }
    inline void addY(float val) { m_vec[1] += val; }
    inline void addZ(float val) { m_vec[2] += val; }
    inline float &operator[](int index) { return m_vec[index]; }
    inline float operator[](int index) const { return m_vec[index]; }
    inline float lengthSquared() { return m_vec[0]*m_vec[0] + m_vec[1]*m_vec[1] + m_vec[2]*m_vec[2]; }
    inline void subtract(const vec3 &v1, const vec3 &v2) { m_vec[0] = v1[0] - v2[0]; m_vec[1] = v1[1] - v2[1]; m_vec[2] = v1[2] - v2[2]; }
private:
    friend std::ostream& operator<<(std::ostream&stream, vec3 vec);
};
