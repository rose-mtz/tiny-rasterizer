#include "MyMath.h"


float determinant(Mat3x3 m)
{
    return (
        m.mat[0][0] * (m.mat[1][1] * m.mat[2][2] - m.mat[1][2] * m.mat[2][1]) -
        m.mat[0][1] * (m.mat[1][0] * m.mat[2][2] - m.mat[1][2] * m.mat[2][0]) + 
        m.mat[0][2] * (m.mat[1][0] * m.mat[2][1] - m.mat[1][1] * m.mat[2][0])
    );
}


Vec3f multiply(Mat3x3 m, Vec3f v)
{
    return Vec3f(
        v.x * m.mat[0][0] + v.y * m.mat[0][1] + v.z * m.mat[0][2],
        v.x * m.mat[1][0] + v.y * m.mat[1][1] + v.z * m.mat[1][2],
        v.x * m.mat[2][0] + v.y * m.mat[2][1] + v.z * m.mat[2][2]
    );
}


// a - b
Vec3f subtract(Vec3f a, Vec3f b)
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


// What the fuck am I doing!