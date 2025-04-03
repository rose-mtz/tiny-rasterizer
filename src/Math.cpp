#include "Math.h"


Vec3f get_triangle_normal(Vec3f a, Vec3f b, Vec3f c)
{
    Vec3f v1 = b - a;
    Vec3f v2 = c - a;

    return (v1 ^ v2).normalize();
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