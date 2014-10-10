#include <cmath>
#include <math/vec3.h>
#include <math/random.h>
namespace CompPhys {
vec3::vec3()
{
    setToZero();
}

vec3::vec3(double x, double y, double z)
{
    set(x,y,z);
}

bool vec3::operator==(vec3 &rhs) {
    return(m_vec[0] == rhs.x() && m_vec[1] == rhs.y() && m_vec[2] == rhs.z());
}

vec3 vec3::operator+(vec3 &rhs) {
    return vec3( m_vec[0] + rhs.x(),
                 m_vec[1] + rhs.y(),
                 m_vec[2] + rhs.z());
}

vec3 vec3::operator-(vec3 &rhs) {
    return vec3( m_vec[0] - rhs.x(),
                 m_vec[1] - rhs.y(),
                 m_vec[2] - rhs.z());
}

vec3 vec3::operator*(vec3 &rhs) {
    return vec3( m_vec[0] * rhs.x(),
                 m_vec[1] * rhs.y(),
                 m_vec[2] * rhs.z());
}

vec3 vec3::operator/(vec3 &rhs) {
    return vec3( m_vec[0] / rhs.x(),
                 m_vec[1] / rhs.y(),
                 m_vec[2] / rhs.z());
}

vec3 vec3::operator+(double scalar) {
    return vec3(m_vec[0] + scalar,
                m_vec[1] + scalar,
                m_vec[2] + scalar);
}

vec3 vec3::operator-(double scalar) {
    return vec3(m_vec[0] - scalar,
                m_vec[1] - scalar,
                m_vec[2] - scalar);
}

vec3 vec3::operator*(double scalar) {
    return vec3(m_vec[0] * scalar,
                m_vec[1] * scalar,
                m_vec[2] * scalar);
}

vec3 vec3::operator/(double scalar) {
    return vec3(m_vec[0] / scalar,
                m_vec[1] / scalar,
                m_vec[2] / scalar);
}

double vec3::dot(vec3 &rhs) {
    return (m_vec[0] * rhs.x() +
            m_vec[1] * rhs.y() +
            m_vec[2] * rhs.z());
}

vec3 vec3::cross(vec3 &rhs) {
    return vec3( m_vec[1] * rhs.z() - m_vec[2] * rhs.y(),
                 m_vec[2] * rhs.x() - m_vec[0] * rhs.z(),
                 m_vec[0] * rhs.y() - m_vec[1] * rhs.x());
}

double vec3::length() {
    return sqrt(lengthSquared());
}

double vec3::lengthSquared() {
    return m_vec[0]*m_vec[0] + m_vec[1]*m_vec[1] + m_vec[2]*m_vec[2];
}

void vec3::normalize() {
    double myLength = length();
    if(myLength > 0) { // Don't divide by zero...
        m_vec[0] /= myLength;
        m_vec[1] /= myLength;
        m_vec[2] /= myLength;
    }
}

void vec3::setToZero()
{
    set(0,0,0);
}

void vec3::randomUniform(double min, double max) {
    m_vec[0] = min + Random::nextDouble()*(max - min);
    m_vec[1] = min + Random::nextDouble()*(max - min);
    m_vec[2] = min + Random::nextDouble()*(max - min);
}

void vec3::set(double x, double y, double z)
{
    m_vec[0] = x;
    m_vec[1] = y;
    m_vec[2] = z;
}

void vec3::randomGaussian(double mean, double standardDeviation) {
    m_vec[0] = Random::nextGaussian(mean, standardDeviation);
    m_vec[1] = Random::nextGaussian(mean, standardDeviation);
    m_vec[2] = Random::nextGaussian(mean, standardDeviation);
}

std::ostream& operator<<(std::ostream &stream, vec3 &vec) {
    return stream << "[" << vec.x() << ", " << vec.y() << ", " << vec.z() << "]";
}
}
