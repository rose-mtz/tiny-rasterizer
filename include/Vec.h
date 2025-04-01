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