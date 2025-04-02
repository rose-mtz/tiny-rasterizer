#include "MyMath.h"
#include <cmath>


float determinant(Mat3x3 m)
{
    return (
        m.mat[0][0] * (m.mat[1][1] * m.mat[2][2] - m.mat[1][2] * m.mat[2][1]) -
        m.mat[0][1] * (m.mat[1][0] * m.mat[2][2] - m.mat[1][2] * m.mat[2][0]) + 
        m.mat[0][2] * (m.mat[1][0] * m.mat[2][1] - m.mat[1][1] * m.mat[2][0])
    );
}


Vec3f operator*(Mat3x3 m, Vec3f v)
{
    return Vec3f(
        v.x * m.mat[0][0] + v.y * m.mat[0][1] + v.z * m.mat[0][2],
        v.x * m.mat[1][0] + v.y * m.mat[1][1] + v.z * m.mat[1][2],
        v.x * m.mat[2][0] + v.y * m.mat[2][1] + v.z * m.mat[2][2]
    );
}


// a - b
Vec3f operator-(Vec3f a, Vec3f b)
{
    return Vec3f(a.x - b.x, a.y - b.y, a.z - b.z);
}


// a x b
Vec3f cross(Vec3f a, Vec3f b)
{
    return Vec3f(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}


float dot(Vec3f a, Vec3f b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


Vec3f operator*(Vec3f v, float s)
{
    return Vec3f(v.x * s, v.y * s, v.z * s);
}


Vec3f operator+(Vec3f a, Vec3f b)
{
    return Vec3f(a.x + b.x, a.y + b.y, a.z + b.z);
}


Vec3f operator/(Vec3f v, float s)
{
    return v * (1.0f / s);
}


float length(Vec3f v)
{
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}


Vec3f normalize(Vec3f v)
{
    return v / length(v);
}


Vec3f get_triangle_normal(Vec3f a, Vec3f b, Vec3f c)
{
    Vec3f v1 = b - a;
    Vec3f v2 = c - a;

    return normalize(cross(v1,v2));
}


Vec3f clampedVec3f(Vec3f v, float min, float max)
{
    return Vec3f(clampedf(v.x, min, max), clampedf(v.y, min, max), clampedf(v.z, min, max));
}


float clampedf(float a, float min, float max)
{
    if (a < min)
    {
        return min;
    }
    else if (a > max)
    {
        return max;
    }
    else
    {
        return a;
    }
}


// What the fuck am I doing!