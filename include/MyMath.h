#pragma once

struct Vec2i
{
    int x, y;

    // Might remove these later
    Vec2i() : x(0), y(0) {}
    Vec2i(int x, int y) : x(x), y(y) {}  
};


struct Vec2f
{
    float x, y;

    // Might remove these later
    Vec2f() : x(0), y(0) {}
    Vec2f(float x, float y) : x(x), y(y) {}  
};


struct Vec3f
{
    float x, y, z;

    // Might remove these later
    Vec3f() : x(0), y(0), z(0) {}
    Vec3f(float a) : x(a), y(a), z(a) {}
    Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
    Vec3f(Vec2i xy, float z) : x(xy.x), y(xy.y), z(z) {}
    Vec3f(Vec2f xy, float z) : x(xy.x), y(xy.y), z(z) {}
};


struct Vec3i
{
    int x, y, z;

    // Might remove these later
    Vec3i() : x(0), y(0), z(0) {}
    Vec3i(int x, int y, int z) : x(x), y(y), z(z) {}  
};


struct Mat3x3
{
    // [row][col], [0][0] is top left
    float mat[3][3];

    Mat3x3() {}
    Mat3x3(Vec3f x, Vec3f y, Vec3f z) 
    {
        mat[0][0] = x.x;
        mat[1][0] = x.y;
        mat[2][0] = x.z;

        mat[0][1] = y.x;
        mat[1][1] = y.y;
        mat[2][1] = y.z;

        mat[0][2] = z.x;
        mat[1][2] = z.y;
        mat[2][2] = z.z;
    }
};

float determinant(Mat3x3 m);
float dot(Vec3f a, Vec3f b);
float length(Vec3f v);
Vec3f cross(Vec3f a, Vec3f b);
Vec3f normalize(Vec3f v);

Vec3f operator*(Mat3x3 m, Vec3f v);
Vec3f operator-(Vec3f a, Vec3f b);
Vec3f operator*(Vec3f v, float s);
Vec3f operator/(Vec3f v, float s);
Vec3f operator+(Vec3f a, Vec3f b);