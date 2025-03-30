#pragma once

struct Vec3f
{
    float x, y, z;

    // Might remove these later
    Vec3f() : x(0), y(0), z(0) {}
    Vec3f(float x, float y, float z) : x(x), y(y), z(z) {}
};


struct Vec3i
{
    int x, y, z;

    // Might remove these later
    Vec3i() : x(0), y(0), z(0) {}
    Vec3i(int x, int y, int z) : x(x), y(y), z(z) {}  
};